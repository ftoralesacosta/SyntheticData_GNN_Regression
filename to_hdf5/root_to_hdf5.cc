#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <H5Cpp.h>
#include <math.h>

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

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
  if (value.GetSetupStatus() < 0) {
    std::cerr << "Error " << value.GetSetupStatus()
      << " setting up reader for " << value.GetBranchName() << '\n';
    return false;
  }
  return true;
}

//Obtain No. Calo hits and MC particles for setting HDF5 dimensions
void find_max_dims(
    char *argv_first[], 
    char *argv_last[], 
    size_t &eventsN_max,
    size_t &hcalNHits_max, 
    size_t &ecalNHits_max, 
    size_t &mcNParticles_max) { 

  for (char **p = argv_first; p != argv_last; p++) {
    TFile *file = TFile::Open(*p);
    TTreeReader events("events", file);

    char hcal_str[] = "HcalEndcapPHitsReco";
    char ecal_str[] = "EcalEndcapPHitsReco";
    char mc_str[]   = "MCParticles";

    TTreeReaderArray<Float_t> hcalE( events, Form("%s.energy",hcal_str));
    TTreeReaderArray<Float_t> ecalE( events, Form("%s.energy",ecal_str));
    TTreeReaderArray<double> mcMass( events, Form("%s.mass",mc_str));

    eventsN_max = events.GetEntries();
    int i = 0;

    fprintf(stderr, "\n%s:%d: %s\n", __func__, __LINE__, "Finding maximum HCal Hits, ECal Hits, and MC Particles from ROOT file");
    while (events.Next()) {
      if (!CheckValue(hcalE)) exit(EXIT_FAILURE);
      /* if (hcalE.GetSize() < 1) continue; */

      size_t hcalNHits = hcalE.GetSize();
      size_t ecalNHits = ecalE.GetSize();
      size_t mcNParticles = mcMass.GetSize();

      hcalNHits_max = std::max(hcalNHits_max, hcalNHits);
      ecalNHits_max = std::max(ecalNHits_max, ecalNHits);
      mcNParticles_max = std::max(mcNParticles_max, mcNParticles);

      if (i%1000 == 0) {
        fprintf(stderr, "%s: %d: event = %i / %u, max HCal Hits= %u, max ECal Hits= %u, max MC Particles = %u\n",
            __func__, __LINE__, i, eventsN_max, hcalNHits_max, ecalNHits_max, mcNParticles_max);
      }

      i++;

    }//while events

    std::cout << std::endl;
    file->Close();
    delete file;

  }//argc

}//find_nhcal_necal_max

void write_data(
    H5::DataSet &hcal_data_set, 
    H5::DataSet &ecal_data_set, 
    H5::DataSet &mc_data_set,
    size_t cal_row_size,
    size_t mc_row_size,
    hsize_t *offset,
    const hsize_t *hcal_dim_extend, 
    const hsize_t *ecal_dim_extend, 
    const hsize_t *mc_dim_extend,
    const UInt_t eventsN_max, 
    const UInt_t hcalNHits_max, 
    const UInt_t ecalNHits_max, 
    const UInt_t mcNParticles_max, 
    const UInt_t block_size,
    char *argv_first[], char *argv_last[]) 
{

  for (char **p = argv_first; p != argv_last; p++) {

    TFile *file = TFile::Open(*p);
    TTreeReader events("events", file);

    char hcal_str[] = "HcalEndcapPHitsReco";
    char ecal_str[] = "EcalEndcapPHitsReco";
    char mc_str[]  =  "MCParticles";

    const char *calo_var[] = {"energy", "position.x", "position.y", "position.z"};
    TTreeReaderArray<Float_t> hcalE( events, Form("%s.energy",hcal_str));
    TTreeReaderArray<Float_t> hcalX( events, Form("%s.position.x",hcal_str));
    TTreeReaderArray<Float_t> hcalY( events, Form("%s.position.y",hcal_str));
    TTreeReaderArray<Float_t> hcalZ( events, Form("%s.position.z",hcal_str));

    TTreeReaderArray<Float_t> ecalE( events, Form("%s.energy",ecal_str));
    TTreeReaderArray<Float_t> ecalX( events, Form("%s.position.x",ecal_str));
    TTreeReaderArray<Float_t> ecalY( events, Form("%s.position.y",ecal_str));
    TTreeReaderArray<Float_t> ecalZ( events, Form("%s.position.z",ecal_str));

    TTreeReaderArray<int> mcPDG( events, Form("%s.PDG",mc_str));
    TTreeReaderArray<int> mcSimulatorStatus( events, Form("%s.simulatorStatus",mc_str));
    TTreeReaderArray<int> mcGeneratorStatus( events, Form("%s.generatorStatus",mc_str));
    TTreeReaderArray<Float_t> mcPX( events, Form("%s.momentum.x",mc_str));
    TTreeReaderArray<Float_t> mcPY( events, Form("%s.momentum.y",mc_str));
    TTreeReaderArray<Float_t> mcPZ( events, Form("%s.momentum.z",mc_str));
    TTreeReaderArray<double> mcMass( events, Form("%s.mass",mc_str));

    //EXAMPLE older Branch Method
    //double hcalE[N_HITS_MAX];
    //events->SetBranchAddress("HcalEndcapPHitsReco.energy", &hcalE);


    std::vector<float> hcal_data(block_size * hcalNHits_max * cal_row_size, NAN);
    std::vector<float> ecal_data(block_size * ecalNHits_max * cal_row_size, NAN);
    std::vector<float> mc_data(block_size * mcNParticles_max * mc_row_size, NAN);

    //Check the max dims are passed correctly

    int i = 0;
    while (events.Next()) {

      int iblock = i % block_size;
      //writing to file is done every [block_size] number of events
      //this variable keeps track of the current increment within a block,
      //as opposed to [i] which is looping through all events

      size_t h_fill = 0; //fill iff non-spikey cell found
      size_t hcalNHits = hcalE.GetSize();

      for (size_t h_hit = 0; h_hit < hcalNHits; h_hit++) 
      {
        if (hcalE[h_hit] > 1e10) continue; //Omit spikey cells
        hcal_data[(iblock*hcalNHits_max + h_fill)*cal_row_size + 0] = hcalE[h_hit] * 1000; //GeV->MeV
        hcal_data[(iblock*hcalNHits_max + h_fill)*cal_row_size + 1] = hcalX[h_hit]; 
        hcal_data[(iblock*hcalNHits_max + h_fill)*cal_row_size + 2] = hcalY[h_hit]; 
        hcal_data[(iblock*hcalNHits_max + h_fill)*cal_row_size + 3] = hcalZ[h_hit]; //note: Z range is 3800-500cm for some reason
        h_fill++;
      }

      size_t e_fill = 0;
      size_t ecalNHits = ecalE.GetSize();
      for (size_t e_hit = 0; e_hit < ecalNHits; e_hit++) 
      {
        if (ecalE[e_hit] > 1e10) continue;
        ecal_data[(iblock*ecalNHits_max + e_fill)*cal_row_size + 0] = ecalE[e_hit] * 1000; //GeV->MeV
        ecal_data[(iblock*ecalNHits_max + e_fill)*cal_row_size + 1] = ecalX[e_hit]; 
        ecal_data[(iblock*ecalNHits_max + e_fill)*cal_row_size + 2] = ecalY[e_hit]; 
        ecal_data[(iblock*ecalNHits_max + e_fill)*cal_row_size + 3] = ecalZ[e_hit]; 
        e_fill++;
      }

      size_t mc_fill = 0;
      size_t mcNParticles = mcMass.GetSize();
      for (size_t particle = 0; particle < mcNParticles; particle++) 
      {
        if (mcGeneratorStatus[particle] != 1) continue;

        //Calculated Values
        float mcPT = hypot(mcPX[particle], mcPY[particle]);
        float mcP = hypot(mcPX[particle], mcPY[particle], mcPZ[particle]);
        float mcTheta = acos(mcPZ[particle]/mcP)*180./M_PI;

        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 0] = mcPDG[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 1] = mcSimulatorStatus[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 2] = mcGeneratorStatus[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 3] = mcPX[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 4] = mcPY[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 5] = mcPZ[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 6] = mcMass[particle]; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 7] = mcPT; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 8] = mcP; 
        mc_data[(iblock*mcNParticles + mc_fill)*mc_row_size + 9] = mcTheta; 
        mc_fill++;
      }// If sim/gun is working, each event should only have one particle with GenStatus = 1

      bool print_hcal = false;
      bool print_mc = false;

      if (iblock == (block_size-1)) {//writes 1 block (2000 events) at a time. Faster/less memory

        if (print_hcal){
          for (size_t h_hit = 0; h_hit < hcalNHits; h_hit++) {
            float E = hcal_data[(iblock*hcalNHits_max + h_hit)*cal_row_size + 0];
            float Z = hcal_data[(iblock*hcalNHits_max + h_hit)*cal_row_size + 3];
            if (std::isnan(E)) break;
            if (h_hit%10 == 0) 
            {
              std::cout <<"HCal Hit # " << h_hit << " / " << hcalNHits << ", E = " << E << std::endl;
              std::cout <<"HCal Hit # " << h_hit << " / " << hcalNHits << ", Z = " << Z << std::endl;
            }
          }
        }

        if (print_mc) {
          for (size_t particle = 0; particle < mcNParticles; particle++) {
            float Mass = mcMass[particle];
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " Generator Status = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 2] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " PX = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 3] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " PY = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 4] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " PZ = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 5] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " Theta = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 9] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles << " Momentum = " << mc_data[(iblock*mcNParticles + particle)*mc_row_size + 8] << std::endl;
            std::cout << std::endl;
            break; //first index should be filled, followed by nans. 
          }           
          std::cout << "Event Number = " << i << std::endl;

        }//print

        if (offset[0] == 0) {
          // Writing the first event. Data spaces is already created with space for one event
          // see file.createDataSet()
          hcal_data_set.write(&hcal_data[0], H5::PredType::NATIVE_FLOAT);
          ecal_data_set.write(&ecal_data[0], H5::PredType::NATIVE_FLOAT);
          mc_data_set.write(&mc_data[0], H5::PredType::NATIVE_FLOAT);

          //Make sure to clear previous arrays
          std::fill(hcal_data.begin(),hcal_data.end(),NAN);
          std::fill(ecal_data.begin(),ecal_data.end(),NAN);
          std::fill(mc_data.begin(),mc_data.end(),NAN);
        }

        else { 
          // Extended-by-1 dimension. First dim is event#
          const hsize_t hcal_dim_extended[RANK] = {
            offset[0] + hcal_dim_extend[0], hcal_dim_extend[1], hcal_dim_extend[2]
          };

          const hsize_t ecal_dim_extended[RANK] = {
            offset[0] + ecal_dim_extend[0], ecal_dim_extend[1], ecal_dim_extend[2]
          };

          const hsize_t mc_dim_extended[RANK] = {
            offset[0] + mc_dim_extend[0], mc_dim_extend[1], mc_dim_extend[2]
          };

          // Extend to the new dimension
          hcal_data_set.extend(hcal_dim_extended);
          ecal_data_set.extend(ecal_dim_extended);
          mc_data_set.extend(mc_dim_extended);

          // Select the hyperslab that only encompasses the
          // difference from extending the data space (i.e. the
          // new event, but offset at the existing event)
          H5::DataSpace hcal_file_space = hcal_data_set.getSpace();
          H5::DataSpace ecal_file_space = ecal_data_set.getSpace();
          H5::DataSpace mc_file_space = mc_data_set.getSpace();

          hcal_file_space.selectHyperslab(
              H5S_SELECT_SET, hcal_dim_extend, offset);

          ecal_file_space.selectHyperslab(
              H5S_SELECT_SET, ecal_dim_extend, offset);

          mc_file_space.selectHyperslab(
              H5S_SELECT_SET, mc_dim_extend, offset);

          // define memory size to fit the extended hyperslab
          H5::DataSpace hcal_memory_space(RANK, hcal_dim_extend, NULL);
          H5::DataSpace ecal_memory_space(RANK, ecal_dim_extend, NULL);
          H5::DataSpace mc_memory_space(RANK, mc_dim_extend, NULL);

          // Write the data from memory space to file space
          hcal_data_set.write(&hcal_data[0], H5::PredType::NATIVE_FLOAT,
              hcal_memory_space, hcal_file_space);

          ecal_data_set.write(&ecal_data[0], H5::PredType::NATIVE_FLOAT,
              ecal_memory_space, ecal_file_space);

          mc_data_set.write(&mc_data[0], H5::PredType::NATIVE_FLOAT,
              mc_memory_space, mc_file_space);

          std::fill(hcal_data.begin(),hcal_data.end(),NAN);
          std::fill(ecal_data.begin(),ecal_data.end(),NAN);
          std::fill(mc_data.begin(),mc_data.end(),NAN);
        }

        offset[0]+=block_size;	 

        std::cout << std::endl;
        fprintf(stderr, "\r%s: %d: %llu / %lld", __func__,__LINE__, offset[0],
            eventsN_max);

      }//BLOCK CHECK
      i++;
    }//Events

    file->Close();
    delete file;
  }//argv files
}//write function

int main(int argc, char *argv[]){
  if (argc < 3) {
    fprintf(stderr, "%s", "Syntax is [command] [root_file] [new hdf5 file name]");
    exit(EXIT_FAILURE);
  }

  static const size_t cal_row_size = 4; //Number of calorimeter hit variables
  static const size_t mc_row_size = 10; //Number of MC truth variables

  size_t eventsN_max = 0;
  size_t hcalNHits_max = 0;
  size_t ecalNHits_max = 0;
  size_t mcNParticles_max = 0;
  size_t block_size = 1000; //affects chunk size, 

  //Constants for Layering HCal
  const double z_offset = 3800.; //[mm]. Fixes some hardcoded setting in ATHENA detector
  const double z_max = 1200.; // actual length of hcal in z [mm]

  find_max_dims(argv + 1, argv + argc - 1, eventsN_max, hcalNHits_max, ecalNHits_max,mcNParticles_max);
  /* eventsN_max = 10000; hcalNHits_max = 1318; ecalNHits_max = 1035; mcNParticles_max = 15; */ 
  //saved for rec_piplus_Energy_0-100GeV.root

  // Access mode H5F_ACC_TRUNC truncates any existing file, while
  // not throwing an exception (unlike H5F_ACC_RDWR)
  std::string file_str = argv[argc-1];
  H5::H5File file( file_str.c_str(), H5F_ACC_TRUNC );
  // The tensor dimension for each new chunk of events

  //The chunking of data can be edited for performance
  hsize_t hcal_dim_extend[RANK] = {block_size, hcalNHits_max, cal_row_size};
  hsize_t ecal_dim_extend[RANK] = {block_size, ecalNHits_max, cal_row_size};
  hsize_t mc_dim_extend[RANK] = {block_size, mcNParticles_max, mc_row_size};

  //Check the hyperslab/extension dimensions are correct
  fprintf(stderr,"\n%s: %d: HDF5 Chunk Cache Size = %u\n",__func__,__LINE__,block_size);
  fprintf(stderr, "%s: %d: hcal_dim_extend = %u %u %u\n", 
      __func__, __LINE__, hcal_dim_extend[0],hcal_dim_extend[1],hcal_dim_extend[2]);
  fprintf(stderr, "%s: %d: ecal_dim_extend = %u %u %u\n", 
      __func__, __LINE__, ecal_dim_extend[0],ecal_dim_extend[1],ecal_dim_extend[2]);
  fprintf(stderr, "%s: %d:  mc_dim_extend  = %u %u  %u\n", 
      __func__, __LINE__, mc_dim_extend[0],mc_dim_extend[1],mc_dim_extend[2]);

  // The maximum tensor dimension, for unlimited number of events
  // a.k.a. the overall dimensions of the dataset
  hsize_t hcal_dim_max[RANK] = {H5S_UNLIMITED, hcalNHits_max, cal_row_size};
  hsize_t ecal_dim_max[RANK] = {H5S_UNLIMITED, ecalNHits_max, cal_row_size};
  hsize_t mc_dim_max[RANK] = {H5S_UNLIMITED, mcNParticles_max, mc_row_size};

  // The extensible HDF5 data space
  H5::DataSpace hcal_data_space(RANK, hcal_dim_extend, hcal_dim_max);
  H5::DataSpace ecal_data_space(RANK, ecal_dim_extend, ecal_dim_max);
  H5::DataSpace mc_data_space(RANK, mc_dim_extend, mc_dim_max);

  // To enable zlib compression (there will be many NANs) and
  // efficient chunking (splitting of the tensor into contingous
  // hyperslabs), a HDF5 property list is needed
  H5::DSetCreatPropList hcal_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList ecal_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList mc_property = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
  // Check for zlib (deflate) availability and enable 
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
      hcal_property.setDeflate(1);
      ecal_property.setDeflate(1);
      mc_property.setDeflate(1);
    }
  }
#endif // HDF5_USE_DEFLATE

  // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
  // the CPU L2 cache size

  hsize_t hcal_dim_chunk[RANK] = {
    hcal_dim_extend[0],
    hcal_dim_extend[1],
    hcal_dim_extend[2]
  };

  hsize_t ecal_dim_chunk[RANK] = {
    ecal_dim_extend[0],
    ecal_dim_extend[1],
    ecal_dim_extend[2]
  };

  hsize_t mc_dim_chunk[RANK] = {
    mc_dim_extend[0],
    mc_dim_extend[1],
    mc_dim_extend[2]
  };

  hcal_property.setChunk(RANK, hcal_dim_chunk);
  ecal_property.setChunk(RANK, ecal_dim_chunk);
  mc_property.setChunk(RANK, mc_dim_chunk);

  // Create the data set, which will have space for the first event (chunk)
  H5::DataSet hcal_data_set =
    file.createDataSet("hcal", H5::PredType::NATIVE_FLOAT,
        hcal_data_space, hcal_property);

  H5::DataSet ecal_data_set =
    file.createDataSet("ecal", H5::PredType::NATIVE_FLOAT,
        ecal_data_space, ecal_property);

  H5::DataSet mc_data_set =
    file.createDataSet("mc", H5::PredType::NATIVE_FLOAT,
        mc_data_space, mc_property);

  hsize_t offset[RANK] = {0, 0, 0};

  write_data(hcal_data_set, ecal_data_set, mc_data_set,
      cal_row_size, mc_row_size,
      offset, hcal_dim_extend, ecal_dim_extend, mc_dim_extend,
      eventsN_max, hcalNHits_max, ecalNHits_max, mcNParticles_max, block_size,argv + 1, argv + argc - 1);
  fprintf(stderr,"\n\n%s: %d: [Complete] \n\n",__func__,__LINE__);

  file.close();

  return EXIT_SUCCESS;
}
