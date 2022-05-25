//This code creates calorimeter 'images' from hdf5 file

// into 3 layers. The two boandaries that define the three
// layers are shifted to various points in the z-direction

#include <iostream>
#include <H5Cpp.h>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace H5;

// Set the CPU L2 cache size to 512 kB
// check /sys/devices/system/cpu/cpu0/cache/index2/size
#ifndef HDF5_DEFAULT_CACHE
#define HDF5_DEFAULT_CACHE (512 * 1024)
#endif

#ifndef HDF5_USE_DEFLATE
#define HDF5_USE_DEFLATE
#endif

#define RANK 3 //Event No., branch, variable index

void get_dims_hdf5(
    hsize_t *dims,
    size_t rank,
    const char *filename, 
    const char *dataset_name)
{
  //Open dataset/dataspace
  H5::H5File h5_file( filename, H5F_ACC_RDONLY );
  H5::DataSet dataset = h5_file.openDataSet(dataset_name);
  H5::DataSpace dataspace = dataset.getSpace();

  //Check number of dimensions, i.e. dataspace rank
  if ( rank != dataspace.getSimpleExtentNdims())
    fprintf(stderr, "\r%s: %d: Dataset RANK mismatch. Check dimensions of dataset.",
        __func__,__LINE__ );

  //Grab Dimensions
  dataspace.getSimpleExtentDims(dims);

  if ( rank == 3)
    fprintf(stderr, "\r%s: %d: %s Dataset Dimensions = [%lu] [%lu] [%lu]\n",
        __func__,__LINE__, dataset_name, dims[0], dims[1], dims[2]);

  else if ( rank == 4)
    fprintf(stderr, "\r%s: %d: %s Dataset Dimensions = [%lu] [%lu] [%lu] [%lu]\n",
        __func__,__LINE__, dataset_name, dims[0], dims[1], dims[2], dims[3]);

  else
    fprintf(stderr, "\r%s: %d: Unsupported Dataset RANK \n",__func__,__LINE__);
}

void add_dataset(
    const char *dset_name,
    size_t rank,
    hsize_t *dims,
    H5::H5File file,
    size_t chunk_size,
    const float FillVal) //FIXME: double check that we don't need to pass by refence. Probably OK as is.
{ //sets Chunking and compression style 

  hsize_t *dim_extend = dims;
  hsize_t *dim_max = dims;
  /* dim_max[0] = H5S_UNLIMITED; */

  H5::DataSpace data_space(rank, dim_extend, dim_max);
  H5::DSetCreatPropList property = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
  if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
    fprintf(stderr, "%s:%d: warning: deflate filter not "
        "available\n", __FILE__, __LINE__);
  }
  else {
    unsigned int filter_info;
    H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
      fprintf(stderr, "%s:%d: warning: deflate filter not "
          "available for encoding\n", __FILE__, __LINE__);
    }
    else {
      property.setDeflate(1);
    }
  }
#endif // HDF5_USE_DEFLATE

  //Set Default Fill Val. NAN for good compression.
  property.setFillValue(H5::PredType::NATIVE_FLOAT, &FillVal);

  //Set Chunk Size [best if set to batch_size from TF]
  hsize_t *dim_chunk = dim_extend; //size_t array to hsize_t array
  dim_chunk[0] = chunk_size;
  property.setChunk(rank, dim_chunk);

  // Create the data set, which will have space for the first event chunk
  H5::DataSet data_set = file.createDataSet(dset_name, 
      H5::PredType::NATIVE_FLOAT, data_space, property);

  return;
}


void create_calo_images(
    const char * calo_h5_name,
    const char * image_h5_name,//FIXME: Change to H5 File as input, not just the name?
    hsize_t* calo_dims,
    hsize_t chunk_size,
    const size_t n_images,
    const size_t n_layers,
    size_t z_offset,
    size_t z_max,
    size_t layer_start = 100,
    size_t layer_max=500,
    size_t z_step = 20)
{
  //Event loop
  //Load ecal and hcal datasets from OLD

  //Open data h5 file
  H5::H5File cal_file(calo_h5_name, H5F_ACC_RDONLY );
  H5::DataSet hcal_dataset = cal_file.openDataSet( "hcal" );
  H5::DataSet ecal_dataset = cal_file.openDataSet( "ecal" );
  H5::DataSpace ecal_dataspace = ecal_dataset.getSpace();
  H5::DataSpace hcal_dataspace = hcal_dataset.getSpace();

  //Open images h5 file
  H5::H5File image_file(image_h5_name, H5F_ACC_RDONLY );
  H5::DataSet image_dataset = image_file.openDataSet( "calo_images" );
  H5::DataSpace image_dataspace = image_dataset.getSpace();

  /* hsize_t* image_dims = get_dims_hdf5(image_h5_name, "calo_images"); */
  hsize_t image_dims[RANK]; 
  get_dims_hdf5(image_dims, RANK, image_h5_name, "calo_images");
  image_dims[0] = chunk_size*n_images;

  size_t N_Events = calo_dims[0];
  size_t N_Hits = calo_dims[2]*2;
  size_t n_variables = image_dims[1];

  float ecal_data[chunk_size][calo_dims[1]][calo_dims[2]] = {NAN};
  float hcal_data[chunk_size][calo_dims[1]][calo_dims[2]] = {NAN};

  size_t rank = hcal_dataspace.getSimpleExtentNdims();
  hsize_t read_dims[rank] = {chunk_size,calo_dims[1],calo_dims[2]};
  hsize_t calo_offset[rank] = {0,0,0}; //offset is incremented in event loop

  //Selection in file
  hcal_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, calo_offset );
  ecal_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, calo_offset );

  //Selection in memeroy (allocate)
  H5::DataSpace calo_memspace(rank, read_dims );
  calo_memspace.selectHyperslab( H5S_SELECT_SET, read_dims, calo_offset );

  //Read from file, into allocated memory
  hcal_dataset.read( hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
  ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);

  hsize_t img_offset[RANK] = {0,0,0};

  std::vector<float>image_data( image_dims[0] * image_dims[1] * image_dims[2], NAN );
  size_t D1 = N_Hits;
  size_t D2 = n_variables;
  size_t D3 = n_images;

  //EVENT LOOP
  for (size_t ievt = 0; ievt < N_Events; ievt++) {
    fprintf(stderr, "\r%s: %d: Processing Event %lu / %lu", __func__,__LINE__,ievt,N_Events );

    size_t ichunk = ievt % chunk_size; //The index within the current event chunk
    size_t i_image = 0;

    //Image Loop 
    for (size_t length_1 = layer_start; length_1 < layer_max; length_1 += z_step) {
      for (size_t length_2 = layer_start; length_2 < layer_max; length_2 += z_step) {

        size_t layer_boundaries[n_layers+1] = 
        {
          z_offset, 
          z_offset+length_1, 
          z_offset+length_1+length_2,
          z_max 
        }; // 4 'edges' defining 3 hcal layers.

        size_t ecal_hit_count = 0;

        //ECal Hit Loop
        for (size_t ihit = 0; ihit < N_Hits; ihit++) {
          if(std::isnan(ecal_data[ichunk][0][ihit])) break;

          size_t x = ihit;
          size_t y = 0; //0-5
          size_t z = i_image;
          size_t t = ichunk; 

          size_t E_index = x + 0 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t X_index = x + 1 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t Y_index = x + 2 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t D_index = x + 3 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t L1_index = x + 4 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t L2_index = x + 5 * D1 + z * D1*D2 + t * D1*D2*D3; 
          //Flat index for '4D' vector [events][images][variable][hit]

          /* fprintf(stderr, "%d: %s: vector index %llu / %llu \n", __LINE__, __func__, index, image_data.size()); */
          size_t ecal_depth = 1;
          image_data[E_index] = hcal_data[ichunk][0][ihit];
          image_data[X_index] = hcal_data[ichunk][1][ihit];
          image_data[Y_index] = hcal_data[ichunk][2][ihit];
          image_data[D_index] = ecal_depth;
          image_data[L1_index] = layer_boundaries[1];
          image_data[L2_index] = layer_boundaries[2];

          ecal_hit_count++;
        }

        //HCal Hit Loop
        for (size_t ihit = 0; ihit < N_Hits; ihit++) {
          if(std::isnan(hcal_data[ichunk][0][ihit])) break;

          size_t x = ihit + ecal_hit_count; //offset by ecal hits
          size_t y = 0; //0-5
          size_t z = i_image;
          size_t t = ichunk; 

          size_t E_index = x + 0 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t X_index = x + 1 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t Y_index = x + 2 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t D_index = x + 3 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t L1_index = x + 4 * D1 + z * D1*D2 + t * D1*D2*D3; 
          size_t L2_index = x + 5 * D1 + z * D1*D2 + t * D1*D2*D3; 
          //Flat index for '4D' vector [events][images][variable][hit]

          /* fprintf(stderr, "%d: %s: vector index %llu / %llu \n", __LINE__, __func__, index, image_data.size()); */
          image_data[E_index] = hcal_data[ichunk][0][ihit];
          image_data[X_index] = hcal_data[ichunk][1][ihit];
          image_data[Y_index] = hcal_data[ichunk][2][ihit];
          image_data[L1_index] = layer_boundaries[1];
          image_data[L2_index] = layer_boundaries[2];

          //Find HCal Depth
          size_t hcal_depth;
          float hcal_z = hcal_data[ichunk][3][ihit];
          for (size_t iz = 0; iz < n_layers; iz++) 
            if ((hcal_z >= layer_boundaries[iz]) && (hcal_z < layer_boundaries[iz+1])) 
              hcal_depth = iz+2;
          /* fprintf(stderr, "%s: %d: [layer][hcalz][boundaries] = [%lu][%1.2f][%lu <-> %lu] \n", */
          /*     __func__, __LINE__,iz+2,hcal_z,layer_boundaries[iz],layer_boundaries[iz+1]); */

          image_data[D_index] = hcal_depth;

          /* hit_count++; */
        }//HCal hit
        i_image++;
      }//layer 1
    }//layer 2

    //Back to Event Loop
    if (ichunk == chunk_size-1){

      /* std::vector<std::vector<float>> image_data; */

      H5::DataSpace img_file_space = image_dataset.getSpace();
      img_file_space.selectHyperslab(H5S_SELECT_SET,image_dims, img_offset);

      // define memory size to fit the extended hyperslab
      H5::DataSpace img_memory_space(RANK,image_dims, NULL);

      image_dataset.write(&image_data[0], H5::PredType::NATIVE_FLOAT,
          img_memory_space, img_file_space);

      std::fill(image_data.begin(),image_data.end(),NAN);
      /* std::fill_n(&image_data[0][0][0][0],chunk_size*n_images*n_variables*N_Hits,NAN); */

      //Read in next chunk of data from calorimeters
      img_offset[0]+=chunk_size*n_images;
      calo_offset[0]+=chunk_size;

      if (img_offset[0] >= N_Events*n_images) break;
      if (calo_offset[0] >= N_Events) break; 

      //Read next set of events
      hcal_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, calo_offset );
      ecal_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, calo_offset );
      hcal_dataset.read(hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
      ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);

    }//event chunk if
    fprintf(stderr, "\r%s: %d: Event %lu / %lu", __func__,__LINE__,ievt,N_Events );
  }//event

  return;
}

void create_truth_data(
    const char *old_hdf5_file,
    const char *new_hdf5_file,
    hsize_t *mc_dims)
{

  //Open MC Dataset
  H5::H5File mc_file(old_hdf5_file, H5F_ACC_RDONLY );
  H5::DataSet mc_dataset = mc_file.openDataSet( "mc" );
  H5::DataSpace mc_dataspace = mc_dataset.getSpace();
  size_t rank = mc_dataspace.getSimpleExtentNdims();

  //Get Dimensions of new 'truth' training dataset
  hsize_t truth_dims[rank]; 
  get_dims_hdf5(truth_dims, rank, new_hdf5_file, "truth");
  size_t chunk_size = truth_dims[0];

  hsize_t read_dims[rank] = {chunk_size,mc_dims[1],mc_dims[2]};
  hsize_t mc_offset[rank] = {0,0,0}; //offset is incremented in event loop

  //Selection in file
  mc_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );

  //Selection in memeroy (allocate)
  H5::DataSpace mc_memspace(rank, read_dims );
  mc_memspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );

  //Read from file, into allocated memory
  float mc_data[chunk_size][mc_dims[1]][mc_dims[2]] = {NAN};
  mc_dataset.read( mc_data, H5::PredType::NATIVE_FLOAT, mc_memspace, mc_dataspace);


  //Open Truth Dataset
  H5::H5File truth_file(new_hdf5_file, H5F_ACC_RDONLY );
  H5::DataSet truth_dataset = truth_file.openDataSet( "truth" );
  H5::DataSpace truth_dataspace = truth_dataset.getSpace();

  std::vector<float> truth_data( chunk_size * truth_dims[1] * truth_dims[2],NAN);
  hsize_t truth_offset[rank] = {0};

  size_t N_Events = mc_dims[0];
  size_t N_Particles = mc_dims[2];
  size_t n_variables = truth_dims[1];

  for (size_t ievt = 0; ievt < N_Events; ievt++) {

    size_t ichunk = ievt % chunk_size; 
    for (size_t particle = 0; particle < N_Particles; particle++) {

      size_t t_P_index = (ichunk*n_variables + 0) * N_Particles + particle;
      size_t t_Theta_index = (ichunk*n_variables + 1) * N_Particles + particle;

      truth_data[t_P_index] = mc_data[ichunk][8][particle]; 
      truth_data[t_Theta_index] = mc_data[ichunk][9][particle];
      //see root_to_hdf5.cc:196 to check indecies 

    }
    if (ichunk == chunk_size-1){

      if (truth_offset[0] == 0) 
        truth_dataset.write(&truth_data[0], H5::PredType::NATIVE_FLOAT);

      else {

        //----------- Write new Truth Data ------------
        // Extended-by-1 dimension. First dim is event#

        const hsize_t truth_dim_extended[RANK] = {
          truth_offset[0] + chunk_size, truth_dims[1], truth_dims[2]};

        truth_dataset.extend(truth_dim_extended);

        H5::DataSpace truth_file_space = truth_dataset.getSpace();
        truth_file_space.selectHyperslab(H5S_SELECT_SET,truth_dims, truth_offset);

        // define memory size to fit the extended hyperslab
        H5::DataSpace truth_memory_space(RANK, truth_dims, NULL);

        // Write the data from memory space to file space
        truth_dataset.write(&truth_data[0], H5::PredType::NATIVE_FLOAT,
            truth_memory_space, truth_file_space);

      }//first chunk else

      truth_offset[0]+=chunk_size;
      std::fill(truth_data.begin(),truth_data.end(),NAN);

      //Read in next chunk of data from calorimeters
      mc_offset[0]+=chunk_size;

      if (mc_offset[0] >= N_Events) break;

      mc_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );
      mc_dataset.read(mc_data, H5::PredType::NATIVE_FLOAT, mc_memspace, mc_dataspace);

    }//chunk check
    fprintf(stderr, "\r%s: %d: Getting Truth Data Event %lu / %lu", __func__,__LINE__,ievt,N_Events );
  }

  return;
}

int main(int argc, char *argv[]){

  if (argc < 3) {
    fprintf(stderr, "%s", "Syntax is [command] [old hdf5 file] [ new layered hdf5 file]");
    exit(EXIT_FAILURE); }

  const char *old_hdf5_file = argv[1];
  const char *new_hdf5_file = argv[2];


  //Grab calorimeter and MC dateset dimensions
  hsize_t calo_dims[RANK];
  hsize_t mc_dims[RANK];
  get_dims_hdf5(calo_dims, RANK, old_hdf5_file, "hcal"); //ecal has same dims
  get_dims_hdf5( mc_dims, RANK, old_hdf5_file, "mc" ); //FIXME: Really only need last element

  //=================New HDF5 File for 'Images'====================//
  //An 'image' is a composition of ecal and hcal data, with various segmenations of the hcal

  H5::H5File image_file( new_hdf5_file, H5F_ACC_TRUNC );

  size_t chunk_events = 100;  //set to tensorflow batch size
  size_t n_truth_vars = 2; //Particle Energy and Theta 
  size_t n_particles_max = mc_dims[2];

  size_t n_img_vars = 6; //[X][Y][E][Depth][Z_Layer1][Z_Layer2] 
  size_t n_images = 400;
  size_t n_hits_max = calo_dims[2];

  size_t n_layers = 3; //No. hcal segmentations
  size_t z_offset = 3800; //[mm]; //FIXME: Remove when generation is updated
  size_t z_max = z_offset + 1200; //HCAL ~ 1.2m
  size_t n_sgmnt_vars = 3;

  hsize_t truth_dims[RANK] = {mc_dims[0], n_truth_vars, mc_dims[2]};
  hsize_t img_dims[RANK] = {calo_dims[0]*n_images, n_img_vars, calo_dims[2]*2};

  const float FillVal = NAN;
  add_dataset("truth", RANK, truth_dims, image_file, chunk_events, FillVal);
  add_dataset("calo_images", RANK, img_dims, image_file, chunk_events, FillVal);

  /* create_truth_data(old_hdf5_file, new_hdf5_file, mc_dims); */
  create_calo_images(old_hdf5_file, new_hdf5_file, calo_dims, chunk_events, n_images, n_layers, z_offset,z_max);  
  //Images: [events X images][variable][n_hits]

  //FIXME: Open H5 file, and pass as argument to functions
  return 0;
}
