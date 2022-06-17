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

#define chunk_events 100 //should match tensorflow batch size

size_t get_rank(H5::DataSet dataset)
{
  H5::DataSpace dataspace = dataset.getSpace();
  return dataspace.getSimpleExtentNdims();
}

void get_dims(H5::DataSet dataset, hsize_t *dims)
{
  H5::DataSpace dataspace = dataset.getSpace();
  dataspace.getSimpleExtentDims(dims);
}

size_t flat_index(
    size_t *dimMax,
    size_t x, 
    size_t y=0, 
    size_t z=0) 
{
  size_t xMax = dimMax[0];
  size_t yMax = dimMax[1];
  size_t zMax = dimMax[2];

  size_t index = x + y*xMax + z*yMax*xMax;

  /* fprintf(stderr, "%s L%d: %llu / %llu, %llu / %llu, %llu / %llu, %llu / %llu\n", */
  /*     __func__,__LINE__, x,xMax, y,yMax, z,zMax, index, xMax*yMax*zMax); */

  return index;
}//gets index for flattend array, up to 3D

size_t get_depth(
    float hcal_z,
    size_t z_offset,
    size_t length_1,
    size_t length_2,
    size_t z_max,
    size_t n_layers=3,
    size_t index_offset=2)
{//index offset for stacking ecal and hcal depths
  
  size_t layer_boundaries[n_layers+1] = 
  {
    z_offset, 
    z_offset+length_1, 
    z_offset+length_1+length_2,
    z_max 
  }; // 4 'edges' defining 3 hcal layers.

  size_t depth;
  for (size_t iz = 0; iz < n_layers; iz++) {
    if ((hcal_z >= layer_boundaries[iz]) && (hcal_z < layer_boundaries[iz+1])) 
    {
      depth = iz+index_offset; 
      //Useful debugg print
      /* fprintf(stderr, "%s %d: iz=%llu depth=%llu, z=%1.2f %llu<->%llu [mm]\n", */
      /*     __func__,__LINE__,iz, depth, hcal_z, layer_boundaries[iz],layer_boundaries[iz+1]); */
    }
  }

  return depth;
}

void add_dataset(
    const char *dset_name,
    hsize_t *dims,
    size_t rank,
    H5::H5File file,
    const float FillVal)
{ //sets Chunking and Compression 

  H5::DataSpace data_space(rank, dims, dims);
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

  //Set Chunk Size
  hsize_t *dim_chunk = dims; 
  dim_chunk[0] = chunk_events; 
  property.setChunk(rank, dim_chunk);
  fprintf(stderr, "%s %d: %s Chunk Size = %llu %llu %llu \n",
      __func__,__LINE__,dset_name,dim_chunk[0],dim_chunk[1],dim_chunk[2]);

  H5::DataSet data_set = file.createDataSet(dset_name, 
      H5::PredType::NATIVE_FLOAT, data_space, property);

  return;
}

void create_calo_images(
    H5::DataSet hcal_dataset,
    H5::DataSet ecal_dataset,
    H5::DataSet image_dataset,
    float Fill_Val,
    const size_t n_images,
    double * img_means,
    double *img_stdevs,
    const size_t n_layers,
    size_t z_offset,
    size_t z_max,
    size_t layer_start = 100,
    size_t layer_max= 500,
    size_t z_step = 20)
{

  //Open DataSpaces
  H5::DataSpace ecal_dataspace = ecal_dataset.getSpace();
  H5::DataSpace hcal_dataspace = hcal_dataset.getSpace();
  H5::DataSpace image_dataspace = image_dataset.getSpace();
  size_t rank = hcal_dataspace.getSimpleExtentNdims();

  hsize_t calo_MaxDims[rank];
  get_dims(hcal_dataset, calo_MaxDims);
  size_t N_Events = calo_MaxDims[0];
  size_t N_calo_vars = calo_MaxDims[1];
  size_t N_calo_hits = calo_MaxDims[2];
  hsize_t calo_chunk_dims[rank] = {chunk_events, N_calo_vars, N_calo_hits};

  hsize_t image_chunk_dims[rank]; 
  get_dims(image_dataset, image_chunk_dims);
  size_t N_img_vars = image_chunk_dims[1];
  size_t N_img_hits = image_chunk_dims[2];
  image_chunk_dims[0] = chunk_events; 


  float ecal_data[chunk_events][N_calo_vars][N_calo_hits] = {NAN};
  float hcal_data[chunk_events][N_calo_vars][N_calo_hits] = {NAN};

  //Selection in file and memory. Read First Chunk
  hsize_t calo_offset[rank] = {0,0,0}; //offset is incremented in event loop
  hcal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  ecal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  H5::DataSpace calo_memspace(rank, calo_chunk_dims );
  calo_memspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  hcal_dataset.read( hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
  ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);

  hsize_t img_offset[rank] = {0,0,0};
  std::vector<float>image_data( chunk_events* N_img_vars * N_img_hits, Fill_Val);
  size_t Img_MaxDims[rank] = {N_img_hits, N_img_vars, chunk_events}; 

  //EVENT LOOP
  /* N_Events = 30; */
  for (size_t ievt = 0; ievt < N_Events; ievt++) {

    size_t ichunk = (ievt % chunk_events); //for reading calo 1 chunk at a time

    //Layer Loop. Each Layer iteration is an image 
    size_t image_count = 0;
    for (size_t length_1 = layer_start; length_1 < layer_max; length_1 += z_step) {
      for (size_t length_2 = layer_start; length_2 < layer_max; length_2 += z_step) {
        size_t i_image = image_count % chunk_events;

        size_t layer_boundaries[n_layers+1] = 
        { z_offset, z_offset+length_1, 
          z_offset+length_1+length_2, z_max };

        //ECal Hit Loop
        size_t ecal_hit_count = 0;
        for (size_t ihit = 0; ihit < N_img_hits; ihit++) {
          if(std::isnan(ecal_data[ichunk][0][ihit])) break;

          size_t ecal_depth = 1;

          size_t E_index  = flat_index(Img_MaxDims, ihit, 0, i_image);
          size_t X_index  = flat_index(Img_MaxDims, ihit, 1, i_image);
          size_t Y_index  = flat_index(Img_MaxDims, ihit, 2, i_image);
          size_t D_index  = flat_index(Img_MaxDims, ihit, 3, i_image);
          size_t L1_index = flat_index(Img_MaxDims, ihit, 4, i_image);
          size_t L2_index = flat_index(Img_MaxDims, ihit, 5, i_image);

          image_data[E_index] = (ecal_data[ichunk][0][ihit] - img_means[0])/img_stdevs[0];
          image_data[X_index] = (ecal_data[ichunk][1][ihit] - img_means[1])/img_stdevs[1];
          image_data[Y_index] = (ecal_data[ichunk][2][ihit] - img_means[2])/img_stdevs[2];
          image_data[D_index] = (ecal_depth - img_means[3])/img_stdevs[3];
          image_data[L1_index] = (layer_boundaries[1] - img_means[4])/img_stdevs[4];
          image_data[L2_index] = (layer_boundaries[2] - img_means[5])/img_stdevs[5];

          ecal_hit_count++;
        }

        //HCal Hit Loop
        for (size_t ihit = 0; ihit < N_img_hits; ihit++) {
          if(std::isnan(hcal_data[ichunk][0][ihit])) break;

          size_t hcal_hit = ihit + ecal_hit_count;
          float  hcal_z = hcal_data[ichunk][3][ihit];
          size_t hcal_depth = get_depth(hcal_z,z_offset,length_1,length_2,z_max);

          size_t E_index  = flat_index(Img_MaxDims, hcal_hit, 0, i_image);
          size_t X_index  = flat_index(Img_MaxDims, hcal_hit, 1, i_image);
          size_t Y_index  = flat_index(Img_MaxDims, hcal_hit, 2, i_image);
          size_t D_index  = flat_index(Img_MaxDims, hcal_hit, 3, i_image);
          size_t L1_index = flat_index(Img_MaxDims, hcal_hit, 4, i_image);
          size_t L2_index = flat_index(Img_MaxDims, hcal_hit, 5, i_image);
          //Flat index for '3D' vector [images][variable][hit]

          image_data[E_index]  = (hcal_data[ichunk][0][ihit] - img_means[0])/img_stdevs[0];
          image_data[X_index]  = (hcal_data[ichunk][1][ihit] - img_means[1])/img_stdevs[1];
          image_data[Y_index]  = (hcal_data[ichunk][2][ihit] - img_means[2])/img_stdevs[2];
          image_data[D_index]  = (hcal_depth - img_means[3])/img_stdevs[3];
          image_data[L1_index] = (layer_boundaries[1] - img_means[4])/img_stdevs[4];
          image_data[L2_index] = (layer_boundaries[2] - img_means[5])/img_stdevs[5];

        }//HCal hit

        //write every 'chunk' number of images
        if (i_image == chunk_events-1){
          // define memory size to fit the extended hyperslab
          H5::DataSpace img_file_space = image_dataset.getSpace();
          img_file_space.selectHyperslab(H5S_SELECT_SET,image_chunk_dims, img_offset);
          H5::DataSpace img_memory_space(rank,image_chunk_dims, NULL);

          image_dataset.write(&image_data[0], H5::PredType::NATIVE_FLOAT,
              img_memory_space, img_file_space);

          std::fill(image_data.begin(),image_data.end(),Fill_Val);

          img_offset[0] += chunk_events; 
          if (img_offset[0] >= N_Events*n_images) break;
          //img dataset is incremented every 'chunk' number of Images
          //but it's total size is N_Events*n_images
        }
        image_count++;//fencepost. Ater set to 0, is added +1
      }//layer 1
    }//layer 2

    //read calo every 'chunk_size' number of events
    if (ichunk == chunk_events-1){

      //Read in next chunk of data from calorimeters
      calo_offset[0]+=chunk_events;
      if (calo_offset[0] >= N_Events) break; 

      //Read next set of events
      hcal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      ecal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      hcal_dataset.read(hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
      ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);

    }//event chunk if

    fprintf(stderr, "\r%s: %d: Getting Image Data Event %lu / %lu", __func__,__LINE__,ievt,N_Events);

  }//event

  fprintf(stderr, "\n%s: %d: Image Data Done \n", __func__,__LINE__ );

  return;
}

void get_mean_stdev(
    double *img_means,
    double *img_stdevs,
    H5::DataSet hcal_dataset,
    H5::DataSet ecal_dataset,
    size_t N_img_vars,
    const size_t n_images,
    const size_t n_layers,
    size_t z_offset,
    size_t z_max,
    size_t layer_start = 100,
    size_t layer_max= 500,
    size_t z_step = 20)
{
  //Open DataSpaces
  H5::DataSpace ecal_dataspace = ecal_dataset.getSpace();
  H5::DataSpace hcal_dataspace = hcal_dataset.getSpace();
  size_t rank = hcal_dataspace.getSimpleExtentNdims();

  hsize_t calo_MaxDims[rank];
  get_dims(hcal_dataset, calo_MaxDims);
  size_t N_Events = calo_MaxDims[0];
  size_t N_calo_vars = calo_MaxDims[1];
  size_t N_calo_hits = calo_MaxDims[2];
  hsize_t calo_chunk_dims[rank] = {chunk_events, N_calo_vars, N_calo_hits};

  float ecal_data[chunk_events][N_calo_vars][N_calo_hits] = {NAN};
  float hcal_data[chunk_events][N_calo_vars][N_calo_hits] = {NAN};

  //Selection in file and memory. Read First Chunk
  hsize_t calo_offset[rank] = {0,0,0}; //offset is incremented in event loop
  hcal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  ecal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  H5::DataSpace calo_memspace(rank, calo_chunk_dims );
  calo_memspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
  hcal_dataset.read( hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
  ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);

  //--------------------Get Mean--------------------
  size_t hit_count = 0;
  for (size_t ievt = 0; ievt < N_Events; ievt++) {
    fprintf(stderr, "\r%s: %d: Calculating Mean and Standard Deviation %lu / %lu",
        __func__,__LINE__,ievt,N_Events );

    size_t ichunk = ievt%chunk_events;

    for (size_t length_1 = layer_start; length_1 < layer_max; length_1 += z_step) {
      for (size_t length_2 = layer_start; length_2 < layer_max; length_2 += z_step) {

        //ECAL DATA
        for (size_t ihit = 0; ihit < N_calo_hits; ihit++) {
          if (isnan(ecal_data[ichunk][0][ihit])) break;
          img_means[0] += ecal_data[ichunk][0][ihit]; //E
          img_means[1] += ecal_data[ichunk][1][ihit]; //X
          img_means[2] += ecal_data[ichunk][2][ihit]; //Y
          img_means[3] += 1; //ecal_depth
          img_means[4] += (z_offset+length_1);
          img_means[5] += (z_offset+length_1+length_2);
          hit_count++;
        }

        //HCAL DATA
        for (size_t ihit = 0; ihit < N_calo_hits; ihit++) {
          if (isnan(hcal_data[ichunk][0][ihit])) break;
          if (isnan(hcal_data[ichunk][3][ihit])) break;

          float hcal_z = hcal_data[ichunk][3][ihit];
          size_t hcal_depth = get_depth(hcal_z,z_offset,length_1,length_2,z_max);

          /* for (size_t ivar = 0; ivar < N_calo_vars; ivar++) { */
          /*   fprintf(stderr, "%s %d: %llu %1.2f\n",__func__,__LINE__,ivar,hcal_data[ichunk][ivar][ihit]); */
          /* } */

          img_means[0] += hcal_data[ichunk][0][ihit]; 
          img_means[1] += hcal_data[ichunk][1][ihit]; 
          img_means[2] += hcal_data[ichunk][2][ihit]; 
          img_means[3] += hcal_depth;
          img_means[4] += (z_offset+length_1);
          img_means[5] += (z_offset+length_1+length_2);
          hit_count++;
        }

      }//L1
    }//L2

    //read in next chunk
    if (ichunk == chunk_events-1) {
      calo_offset[0] += chunk_events; 
      if (calo_offset[0] >= N_Events) break; 
      hcal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      ecal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      hcal_dataset.read(hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
      ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);
    }

  }//event loop
  fprintf(stderr, "\n");

  for (size_t ivar = 0; ivar < N_img_vars; ivar++)
    img_means[ivar] = img_means[ivar] / hit_count; //sum / N

  //--------------------Get Standard Deviation--------------------
  for (size_t ievt = 0; ievt < N_Events; ievt++) {
    size_t ichunk = ievt%chunk_events;
    //ecal
    for (size_t ihit = 0; ihit < N_calo_hits; ihit++) {
      if (isnan(ecal_data[ichunk][0][ihit])) break;
      for (size_t ivar = 0; ivar < N_calo_vars; ivar++)
        img_stdevs[ivar] += pow(ecal_data[ichunk][ivar][ihit] - img_means[ivar], 2);
      for (size_t length_1 = layer_start; length_1 < layer_max; length_1 += z_step) {
        for (size_t length_2 = layer_start; length_2 < layer_max; length_2 += z_step) {
          img_stdevs[3] += pow(1 - img_means[3], 2); //ecal_depth
          img_stdevs[4] += pow((z_offset+length_1) - img_means[4], 2);
          img_stdevs[5] += pow((z_offset+length_1+length_2) - img_means[5], 2);
        }
      }
    }
    //hcal
    for (size_t ihit = 0; ihit < N_calo_hits; ihit++) {
      if (isnan(hcal_data[ichunk][0][ihit])) break;
      for (size_t ivar = 0; ivar < N_calo_vars; ivar++)
        img_stdevs[ivar] += pow(hcal_data[ichunk][ivar][ihit] - img_means[ivar], 2);

      for (size_t length_1 = layer_start; length_1 < layer_max; length_1 += z_step) {
        for (size_t length_2 = layer_start; length_2 < layer_max; length_2 += z_step) {

          size_t layer_boundaries[n_layers+1] = 
          { z_offset, z_offset+length_1, 
            z_offset+length_1+length_2, z_max };

          float hcal_z = hcal_data[ichunk][3][ihit];
          for (size_t iz = 0; iz < n_layers; iz++) 
            if ((hcal_z >= layer_boundaries[iz]) && (hcal_z < layer_boundaries[iz+1])) 
              img_stdevs[3] += pow((iz+2) - img_means[3], 2);
          //FIXME: Use the depth function for better code standards

          img_stdevs[4] += pow(layer_boundaries[1] - img_means[4], 2);
          img_stdevs[5] += pow(layer_boundaries[2] - img_means[5], 2);
        }
      }
    }
    //Read in next chunk
    if (ichunk == chunk_events-1) {
      calo_offset[0] += chunk_events; 
      if (calo_offset[0] >= N_Events) break; 
      hcal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      ecal_dataspace.selectHyperslab( H5S_SELECT_SET, calo_chunk_dims, calo_offset );
      hcal_dataset.read(hcal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, hcal_dataspace);
      ecal_dataset.read( ecal_data, H5::PredType::NATIVE_FLOAT, calo_memspace, ecal_dataspace);
    }
  }

  for (size_t ivar = 0; ivar < N_img_vars; ivar++){
    if (ivar <= 2)
      img_stdevs[ivar] = std::sqrt(img_stdevs[ivar] / hit_count);
    else 
      img_stdevs[ivar] = std::sqrt(img_stdevs[ivar] / hit_count / n_images);
    fprintf(stderr, "%s: %d: Variable %llu Mean = %1.2f StDev = %1.2f\n", __func__, __LINE__, ivar, img_means[ivar],img_stdevs[ivar]);
  }

}

void create_truth_data(
    H5::DataSet mc_dataset,
    H5::DataSet truth_dataset,
    float Fill_Val,
    size_t n_images)
{

  //Open DataSpaces
  H5::DataSpace mc_dataspace = mc_dataset.getSpace();
  H5::DataSpace truth_dataspace = truth_dataset.getSpace();
  size_t rank = mc_dataspace.getSimpleExtentNdims();

  //Get Dataset Dimensions
  hsize_t mc_dims[rank];
  get_dims(mc_dataset, mc_dims);

  hsize_t truth_dims[rank]; 
  get_dims(truth_dataset,truth_dims);
  truth_dims[0] = chunk_events;

  //Selection in memeroy
  hsize_t read_dims[rank] = {chunk_events,mc_dims[1],mc_dims[2]};
  hsize_t mc_offset[rank] = {0,0,0}; //offset is incremented in event loop
  mc_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );
  H5::DataSpace mc_memspace(rank, read_dims);
  mc_memspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );
  //Read from file, into allocated memory
  float mc_data[chunk_events][mc_dims[1]][mc_dims[2]] = {NAN};
  mc_dataset.read( mc_data, H5::PredType::NATIVE_FLOAT,mc_dataspace, mc_dataspace);

  size_t N_Events = mc_dims[0];
  size_t n_variables = truth_dims[1];
  size_t N_Particles = mc_dims[2];
  size_t MaxDims[rank] = {N_Particles,n_variables,n_images};

  std::vector<float> truth_data(chunk_events * n_variables * N_Particles, Fill_Val);
  hsize_t truth_offset[rank] = {0};

  /* N_Events = 30; */
  for (size_t ievt = 0; ievt < N_Events; ievt++) {

    size_t ichunk = ievt % chunk_events; 

    for (size_t i_image = 0; i_image < n_images; i_image++) {

      size_t img_chunk = i_image % chunk_events;

      for (size_t particle = 0; particle < N_Particles; particle++) {

        if (std::isnan(mc_data[ichunk][8][particle])) continue;
        size_t P_index = flat_index(MaxDims, particle, 0, img_chunk);
        size_t Theta_index = flat_index(MaxDims, particle, 1, img_chunk);
        truth_data[P_index] = mc_data[ichunk][8][particle]; 
        truth_data[Theta_index] = mc_data[ichunk][9][particle];
        //see root_to_hdf5.cc:196 to check indecies 

      }

      //----------- Write new Truth Data ------------
      if (img_chunk == chunk_events-1){
        H5::DataSpace truth_file_space = truth_dataset.getSpace();
        truth_file_space.selectHyperslab(H5S_SELECT_SET,truth_dims, truth_offset);
        H5::DataSpace truth_memory_space(rank, truth_dims, NULL);

        truth_dataset.write(&truth_data[0], H5::PredType::NATIVE_FLOAT,
            truth_memory_space, truth_file_space);

        truth_offset[0]+=chunk_events;
        if (truth_offset[0] >= N_Events*n_images) break;

      }
    }//image loop

    if (ichunk == chunk_events-1){
      // Extended-by-1 dimension. First dim is event#
      mc_offset[0]+=chunk_events;
      std::fill(truth_data.begin(),truth_data.end(),Fill_Val);

      if (mc_offset[0] >= N_Events) break;
      mc_dataspace.selectHyperslab( H5S_SELECT_SET, read_dims, mc_offset );
      mc_dataset.read(mc_data, H5::PredType::NATIVE_FLOAT,mc_memspace, mc_dataspace);

    }//chunk check
    fprintf(stderr, "\r%s: %d: Getting Truth Data Event %lu / %lu", __func__,__LINE__,ievt,N_Events );
  }
  fprintf(stderr, "\n%s: %d: Truth Data Done \n", __func__,__LINE__ );

  return;
}

int main(int argc, char *argv[]){

  if (argc < 3) {
    fprintf(stderr, "%s", "Syntax is [command] [old hdf5 file] [ name of new layered hdf5 file]");
    exit(EXIT_FAILURE); }

  const char *old_hdf5_file = argv[1];
  const char *new_hdf5_file = argv[2];

  //Grab Calorimeter and MC Dataspaces
  H5::H5File calo_file(old_hdf5_file, H5F_ACC_RDONLY );
  H5::DataSet hcal_dataset = calo_file.openDataSet( "hcal" );
  H5::DataSet ecal_dataset = calo_file.openDataSet( "ecal" );
  H5::DataSet  mc_dataset  = calo_file.openDataSet(  "mc"  );

  //Grab calorimeter and MC dateset dimensions
  size_t rank = get_rank(hcal_dataset);
  hsize_t calo_dims[rank];
  hsize_t mc_dims[rank];
  get_dims(hcal_dataset, calo_dims);
  get_dims(mc_dataset, mc_dims);

  //==================================New HDF5 File for 'Images'==================================//
  //An 'image' is a composition of ecal and hcal data, with various segmentations of the hcal
  H5::H5File image_file( new_hdf5_file, H5F_ACC_TRUNC );

  //Parameters for Images/layers, and dataset for training  
  size_t N_Events = calo_dims[0];
  /* N_Events = 30; */

  size_t n_truth_vars = 2; //Particle Energy and Theta 
  size_t n_particles_max = mc_dims[2];

  size_t n_img_vars = 6; //[X][Y][E][Depth][Z_Layer1][Z_Layer2] 
  size_t n_images = 400;
  size_t n_hits_max = calo_dims[2];

  size_t n_layers = 3; //No. hcal segmentations
  size_t z_offset = 3800; //[mm]; //FIXME: Remove when generation is updated
  size_t z_max = z_offset + 1200; //HCAL ~ 1.2m
  size_t n_sgmnt_vars = 3;

  //Create new Image and Truth datasets
  hsize_t truth_dims[rank] = {N_Events*n_images, n_truth_vars, n_particles_max};
  hsize_t img_dims[rank] = {N_Events*n_images, n_img_vars, calo_dims[2]*2}; //*2: hcal+ecal hits

  const float FillVal = 0; //0 instead of NAN for easier regression
  add_dataset("truth", truth_dims, rank, image_file, FillVal);
  add_dataset("calo_images", img_dims, rank, image_file, FillVal);

  H5::DataSet image_dataset = image_file.openDataSet( "calo_images" );
  H5::DataSet truth_dataset = image_file.openDataSet( "truth" );

  //Get Mean and Stdev for each variable
  double img_means[n_img_vars] = {0};
  double img_stdevs[n_img_vars] = {0};

  get_mean_stdev(img_means, img_stdevs, hcal_dataset, 
      ecal_dataset, n_img_vars, n_images, n_layers, z_offset, z_max);
  //image norm done here for performance
  //truth norm done in jupyter for convenience

  create_calo_images(hcal_dataset, ecal_dataset, 
      image_dataset, FillVal, n_images, img_means, img_stdevs, n_layers, z_offset,z_max);  

  create_truth_data(mc_dataset, truth_dataset, FillVal, n_images);

  return 0;
}
