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

#define MAX_SUBSYS 10 // max number of subsystems

bool verbose = true ;

int n_subsystems ;
char subsystem_prefixes[MAX_SUBSYS][100]  ;
char subsystem_short_names[MAX_SUBSYS][10] ;

//-----------------------------------------------------------------------------------------------------

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
  if (value.GetSetupStatus() < 0) {
    std::cerr << "Error " << value.GetSetupStatus()
      << " setting up reader for " << value.GetBranchName() << '\n';
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------------------------------

//Obtain No. Calo hits and MC particles for setting HDF5 dimensions
void find_max_dims(
    char *argv_first[],
    char *argv_last[],
    size_t &eventsN_max,
    size_t &calo_NHits_max,
    size_t &mcNParticles_max) {


  for (char **p = argv_first; p != argv_last; p++) {

    TFile *file = TFile::Open(*p);

    printf("  find_max_dims : Opening %s\n", *p ) ;

    TTreeReader events("events", file);



    eventsN_max = events.GetEntries();
    //eventsN_max = 20000; // FOR TESTING






    events.Restart() ;
    TTreeReaderArray<double> mcMass( events, Form("MCParticles.mass"));
    TTreeReaderArray<Float_t>* subsysE[MAX_SUBSYS] ;
    for ( int si=0; si<n_subsystems; si++ ) {
          subsysE[si] = new TTreeReaderArray<Float_t>( events, Form("%s.energy", subsystem_prefixes[si]));
    }

    {
       int i = 0;
       while (events.Next()) {

          if (i >= eventsN_max) break; // FOR TESTING

          size_t mcNParticles = mcMass.GetSize();

          mcNParticles_max = std::max(mcNParticles_max, mcNParticles);

          for ( int si=0; si<n_subsystems; si++ ) {
               size_t subsysNHits = subsysE[si] -> GetSize();
               calo_NHits_max = std::max(calo_NHits_max, subsysNHits);
          }

          if (i%1000 == 0) {
            fprintf(stderr, "%s: %d: event = %i / %u, max Calo Hits= %u, max MC Particles = %u\n",
                __func__, __LINE__, i, eventsN_max, calo_NHits_max, mcNParticles_max);
          }

          i++;

       }//while events
    }
    for ( int si=0; si<n_subsystems; si++ ) {
          delete subsysE[si] ;
    }




    std::cout << std::endl;
    file->Close();
    delete file;



    printf("\n\n  find_max_dims :   eventsN_max = %u,  calo_NHits_max = %u,  mcNParticles_max = %u\n\n",
       eventsN_max, calo_NHits_max, mcNParticles_max ) ;

  }//argc

}//find_max_dims

//-----------------------------------------------------------------------------------------------------

//void write_data(
//    H5::DataSet &hcal_data_set, 
//    H5::DataSet &ecal_data_set, 
//    H5::DataSet &mc_data_set,
//    size_t cal_row_size,
//    size_t mc_row_size,
//    hsize_t *offset,
//    const hsize_t *hcal_dim_extend, 
//    const hsize_t *ecal_dim_extend, 
//    const hsize_t *mc_dim_extend,
//    const UInt_t eventsN_max, 
//    const UInt_t calo_NHits_max, 
//    const UInt_t mcNParticles_max, 
//    const UInt_t block_size,
//    char *argv_first[], char *argv_last[]) 


void write_data(
    H5::DataSet* calo_data_set[],
    H5::DataSet &mc_data_set,
    size_t cal_row_size,
    size_t mc_row_size,
    hsize_t *offset,
    const hsize_t calo_dim_extend[][RANK],
    const hsize_t *mc_dim_extend,
    const UInt_t eventsN_max,
    const UInt_t calo_NHits_max,
    const UInt_t mcNParticles_max,
    const UInt_t block_size,
    char *argv_first[], char *argv_last[]) 
{

  printf("\n\n ========== Beginning of write_data\n\n") ;

  for (char **p = argv_first; p != argv_last; p++) {

    TFile *file = TFile::Open(*p);
    TTreeReader events("events", file);

//  char hcal_str[] = "HcalEndcapPHitsReco";
//  char ecal_str[] = "EcalEndcapPHitsReco";

    char mc_str[]  =  "MCParticles";

//  const char *calo_var[] = {"energy", "position.x", "position.y", "position.z"};
//  TTreeReaderArray<Float_t> hcalE( events, Form("%s.energy",hcal_str));
//  TTreeReaderArray<Float_t> hcalT( events, Form("%s.time",hcal_str));
//  TTreeReaderArray<Float_t> hcalX( events, Form("%s.position.x",hcal_str));
//  TTreeReaderArray<Float_t> hcalY( events, Form("%s.position.y",hcal_str));
//  TTreeReaderArray<Float_t> hcalZ( events, Form("%s.position.z",hcal_str));

//  TTreeReaderArray<Float_t> ecalE( events, Form("%s.energy",ecal_str));
//  TTreeReaderArray<Float_t> ecalT( events, Form("%s.time",ecal_str));
//  TTreeReaderArray<Float_t> ecalX( events, Form("%s.position.x",ecal_str));
//  TTreeReaderArray<Float_t> ecalY( events, Form("%s.position.y",ecal_str));
//  TTreeReaderArray<Float_t> ecalZ( events, Form("%s.position.z",ecal_str));


    TTreeReaderArray<Float_t>* subsysE[MAX_SUBSYS] ;
    TTreeReaderArray<Float_t>* subsysT[MAX_SUBSYS] ;
    TTreeReaderArray<Float_t>* subsysX[MAX_SUBSYS] ;
    TTreeReaderArray<Float_t>* subsysY[MAX_SUBSYS] ;
    TTreeReaderArray<Float_t>* subsysZ[MAX_SUBSYS] ;
    for ( int si=0; si<n_subsystems; si++ ) {
          subsysE[si] = new TTreeReaderArray<Float_t>( events, Form("%s.energy", subsystem_prefixes[si]));
          subsysT[si] = new TTreeReaderArray<Float_t>( events, Form("%s.time", subsystem_prefixes[si]));
          subsysX[si] = new TTreeReaderArray<Float_t>( events, Form("%s.position.x", subsystem_prefixes[si]));
          subsysY[si] = new TTreeReaderArray<Float_t>( events, Form("%s.position.y", subsystem_prefixes[si]));
          subsysZ[si] = new TTreeReaderArray<Float_t>( events, Form("%s.position.z", subsystem_prefixes[si]));
    } // si




    TTreeReaderArray<int> mcPDG( events, Form("%s.PDG",mc_str));
    TTreeReaderArray<int> mcSimulatorStatus( events, Form("%s.simulatorStatus",mc_str));
    TTreeReaderArray<int> mcGeneratorStatus( events, Form("%s.generatorStatus",mc_str));
    TTreeReaderArray<Float_t> mcPX( events, Form("%s.momentum.x",mc_str));
    TTreeReaderArray<Float_t> mcPY( events, Form("%s.momentum.y",mc_str));
    TTreeReaderArray<Float_t> mcPZ( events, Form("%s.momentum.z",mc_str));
    TTreeReaderArray<double> mcMass( events, Form("%s.mass",mc_str));


//  std::vector<float> hcal_data(block_size * cal_row_size * calo_NHits_max, NAN);
//  std::vector<float> ecal_data(block_size * cal_row_size * calo_NHits_max, NAN);

    std::vector<float>* subsys_data[MAX_SUBSYS] ;
    for ( int si=0; si<n_subsystems; si++ ) {
       subsys_data[si] = new std::vector<float>( block_size * cal_row_size * calo_NHits_max, NAN );
    } // si


    std::vector<float> mc_data(block_size *  mc_row_size * mcNParticles_max, NAN);

    //Check the max dims are passed correctly

    int i = 0;
    while (events.Next()) {

      if (i >= eventsN_max) break;
      int iblock = i % block_size;
      //writing to file is done every [block_size] number of events
      //this variable keeps track of the current increment within a block,
      //as opposed to [i] which is looping through all events. Note 'i' is 
      //incremented only after an event passes selection (hit requirements)


//    size_t h_fill = 0; //fill iff non-spikey cell found
//    size_t hcalNHits = hcalE.GetSize();


      size_t subsys_fill[MAX_SUBSYS] ;
      for ( int si=0; si<n_subsystems; si++ ) { subsys_fill[si] = 0 ; }

      size_t subsysNHits[MAX_SUBSYS] ;
      for ( int si=0; si<n_subsystems; si++ ) { subsysNHits[si] = subsysE[si] -> GetSize() ; }



//    for (size_t h_hit = 0; h_hit < hcalNHits; h_hit++) 
//    {
//      if (hcalE[h_hit] > 1e10) continue; //Omit spikey cells
//      if (hcalE[h_hit] <= 0.00006 ) continue; //Omit Empty Cells and MIPS
//      if (hcalT[h_hit] > 200) continue; //cut long tails in time (realistic for EIC)
//      
//      size_t E_index = iblock*cal_row_size*calo_NHits_max + 0*calo_NHits_max + h_fill;
//      size_t X_index = iblock*cal_row_size*calo_NHits_max + 1*calo_NHits_max + h_fill;
//      size_t Y_index = iblock*cal_row_size*calo_NHits_max + 2*calo_NHits_max + h_fill;
//      size_t Z_index = iblock*cal_row_size*calo_NHits_max + 3*calo_NHits_max + h_fill;
//      //Index for flattened 3D vector
//      
//      hcal_data[E_index] = hcalE[h_hit] *1000;
//      hcal_data[X_index] = hcalX[h_hit]; 
//      hcal_data[Y_index] = hcalY[h_hit]; 
//      hcal_data[Z_index] = hcalZ[h_hit]; 

//      //OLD WAY
//      /* hcal_data[(iblock*cal_row_size + 1)*calo_NHits_max + h_fill] = hcalX[h_hit]; */ 

//      h_fill++;
//    }
//    if (hcalNHits == 0) continue; //skip events with no hcal hits
//    if (h_fill == 0) continue; //skip events hits passing selection

      

//    size_t e_fill = 0;
//    size_t ecalNHits = ecalE.GetSize();

//    for (size_t e_hit = 0; e_hit < ecalNHits; e_hit++) 
//    {
//      if (ecalE[e_hit] > 1e10) continue;
//      if (ecalE[e_hit] <= 0.00006 ) continue; 
//      if (ecalT[e_hit] > 200 ) continue; 

//      size_t E_index = iblock*cal_row_size*calo_NHits_max + 0*calo_NHits_max + e_fill;
//      size_t X_index = iblock*cal_row_size*calo_NHits_max + 1*calo_NHits_max + e_fill;
//      size_t Y_index = iblock*cal_row_size*calo_NHits_max + 2*calo_NHits_max + e_fill;
//      size_t Z_index = iblock*cal_row_size*calo_NHits_max + 3*calo_NHits_max + e_fill;
//      //Index for flattened 3D vector
//      
//      ecal_data[E_index] = ecalE[e_hit] *1000;
//      ecal_data[X_index] = ecalX[e_hit]; 
//      ecal_data[Y_index] = ecalY[e_hit]; 
//      ecal_data[Z_index] = ecalZ[e_hit]; 

//      e_fill++;
//    }
//    if (ecalNHits == 0) continue; //skip events with no hcal hits
//    if (e_fill == 0) continue; //skip events hits passing selection




      for ( int si=0; si<n_subsystems; si++ ) {
         for (size_t h_hit = 0; h_hit < subsysNHits[si]; h_hit++) {
            if ( (*(subsysE[si]))[h_hit] > 1e10) continue ;
            if ( (*(subsysE[si]))[h_hit] <= 0.00006 ) continue ;
            if ( (*(subsysT[si]))[h_hit] > 200 ) continue ;
            size_t E_index = iblock*cal_row_size*calo_NHits_max + 0*calo_NHits_max + subsys_fill[si];
            size_t X_index = iblock*cal_row_size*calo_NHits_max + 1*calo_NHits_max + subsys_fill[si];
            size_t Y_index = iblock*cal_row_size*calo_NHits_max + 2*calo_NHits_max + subsys_fill[si];
            size_t Z_index = iblock*cal_row_size*calo_NHits_max + 3*calo_NHits_max + subsys_fill[si];
            //Index for flattened 3D vector
            (*(subsys_data[si]))[E_index] = (*(subsysE[si]))[h_hit] * 1000 ;
            (*(subsys_data[si]))[X_index] = (*(subsysX[si]))[h_hit] ;
            (*(subsys_data[si]))[Y_index] = (*(subsysY[si]))[h_hit] ;
            (*(subsys_data[si]))[Z_index] = (*(subsysZ[si]))[h_hit] ;
            subsys_fill[si] ++ ;
         } // h_hit
         if ( subsysNHits[si] == 0 ) continue ;
         if ( subsys_fill[si] == 0 ) continue ;
      } // si


      size_t mc_fill = 0;
      size_t mcNParticles = mcMass.GetSize();
      for (size_t particle = 0; particle < mcNParticles; particle++) 
      {
        if (mcGeneratorStatus[particle] != 1) continue;

        //Calculated Values
        float mcPT = hypot(mcPX[particle], mcPY[particle]);
        float mcP = hypot(mcPX[particle], mcPY[particle], mcPZ[particle]);
        float mcTheta = acos(mcPZ[particle]/mcP)*180./M_PI;

        mc_data[(iblock*mc_row_size + 0)*mcNParticles_max + mc_fill] = mcPDG[particle]; 
        mc_data[(iblock*mc_row_size + 1)*mcNParticles_max + mc_fill] = mcSimulatorStatus[particle]; 
        mc_data[(iblock*mc_row_size + 2)*mcNParticles_max + mc_fill] = mcGeneratorStatus[particle]; 
        mc_data[(iblock*mc_row_size + 3)*mcNParticles_max + mc_fill] = mcPX[particle]; 
        mc_data[(iblock*mc_row_size + 4)*mcNParticles_max + mc_fill] = mcPY[particle]; 
        mc_data[(iblock*mc_row_size + 5)*mcNParticles_max + mc_fill] = mcPZ[particle]; 
        mc_data[(iblock*mc_row_size + 6)*mcNParticles_max + mc_fill] = mcMass[particle]; 
        mc_data[(iblock*mc_row_size + 7)*mcNParticles_max + mc_fill] = mcPT; 
        mc_data[(iblock*mc_row_size + 8)*mcNParticles_max + mc_fill] = mcP; 
        mc_data[(iblock*mc_row_size + 9)*mcNParticles_max + mc_fill] = mcTheta; 

        /* std::cout << " " << std::endl; */
        /* for (size_t ivar = 0; ivar < mc_row_size; ivar++) { */
        /*   std::cout << "index = " << (iblock*mc_row_size + ivar)*mcNParticles_max + mc_fill<< " / " << mc_data.size() << std::endl; */
        /* } */
        /* std::cout << " " << std::endl; */

        mc_fill++;
      }// If sim/gun is working, each event should only have one particle with GenStatus = 1
      if (mc_fill == 0) continue;

      bool print_cal = false;
      bool print_mc = false;
      //bool print_cal = true;
      //bool print_mc = true;



      if (iblock == (block_size-1)) 
      {
        //writes 1 block (100 events) at a time. Faster/less memory

        if (print_cal){
          for ( int si=0; si<n_subsystems; si++ ) {
             for (size_t h_hit = 0; h_hit < subsysNHits[si]; h_hit++) {
               //float E = hcal_data[(iblock*cal_row_size + 0)*calo_NHits_max + h_hit];
               //float Z = hcal_data[(iblock*cal_row_size + 3)*calo_NHits_max + h_hit];
               float E = (*(subsys_data[si]))[(iblock*cal_row_size + 0)*calo_NHits_max + h_hit];
               float Z = (*(subsys_data[si]))[(iblock*cal_row_size + 3)*calo_NHits_max + h_hit];
               if (std::isnan(E)) break;
               if (h_hit%10 == 0)
               {
                 //std::cout <<"HCal Hit # " << h_hit << " / " << hcalNHits << ", E = " << E << std::endl;
                 //std::cout <<"HCal Hit # " << h_hit << " / " << hcalNHits << ", Z = " << Z << std::endl;
                 printf("   %s Hit # %5d / %5d, E = %10.6f\n", subsystem_prefixes[si], h_hit, subsysNHits[si], E ) ;
                 printf("   %s Hit # %5d / %5d, Z = %10.6f\n", subsystem_prefixes[si], h_hit, subsysNHits[si], Z ) ;
               }
             } // h_hit
          } // si
        } // print_cal?

        if (print_mc) {
          for (size_t particle = 0; particle < mcNParticles_max; particle++) {
            float Mass = mcMass[particle];
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " GenStatus = "<<mc_data[(iblock*mc_row_size+2)*mcNParticles_max+particle]<<std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " PX = " << mc_data[(iblock*mc_row_size + 3)*mcNParticles_max +particle] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " PY = " << mc_data[(iblock*mc_row_size + 4)*mcNParticles_max +particle] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " PZ = " << mc_data[(iblock*mc_row_size + 5)*mcNParticles_max +particle] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " Theta = "<< mc_data[(iblock*mc_row_size+9)*mcNParticles_max +particle] << std::endl;
            std::cout << "MC Particle # " << particle << " / " << mcNParticles_max << " Momentum = "<< mc_data[(iblock*mc_row_size+8)*mcNParticles_max+particle]<<std::endl;
            std::cout << std::endl;
            break; //first index should be filled, followed by nans. rm this if want to see all particles
          }           
          std::cout << "Event Number = " << i << std::endl;

        }//print

        if (offset[0] == 0) {
          // Writing the first event. Data spaces is already created with space for one event
          // see file.createDataSet()
          ////////hcal_data_set.write(&hcal_data[0], H5::PredType::NATIVE_FLOAT);
          ////////ecal_data_set.write(&ecal_data[0], H5::PredType::NATIVE_FLOAT);
          for ( size_t si=0; si<n_subsystems; si++ ) {
             calo_data_set[si] -> write( &(*(subsys_data[si]))[0], H5::PredType::NATIVE_FLOAT );
          } // si
          mc_data_set.write(&mc_data[0], H5::PredType::NATIVE_FLOAT);

          //Make sure to clear previous arrays
          ////////////std::fill(hcal_data.begin(),hcal_data.end(),NAN);
          ////////////std::fill(ecal_data.begin(),ecal_data.end(),NAN);
          for ( size_t si=0; si<n_subsystems; si++ ) {
             std::fill( (*(subsys_data[si])).begin(), (*(subsys_data[si])).end(),NAN);
          }
          std::fill(mc_data.begin(),mc_data.end(),NAN);
        }

        else {

//        // Extended-by-1 dimension. First dim is event#
//        const hsize_t hcal_dim_extended[RANK] = 
//        {
//          offset[0] + hcal_dim_extend[0], hcal_dim_extend[1], hcal_dim_extend[2]
//        };

//        const hsize_t ecal_dim_extended[RANK] = 
//        {
//          offset[0] + ecal_dim_extend[0], ecal_dim_extend[1], ecal_dim_extend[2]
//        };


          hsize_t subsys_dim_extended[MAX_SUBSYS][RANK] ;

          for ( size_t si=0; si<n_subsystems; si++ ) {
             subsys_dim_extended[si][0] = offset[0] + calo_dim_extend[si][0] ;
             subsys_dim_extended[si][1] = calo_dim_extend[si][1] ;
             subsys_dim_extended[si][2] = calo_dim_extend[si][2] ;
          } // si


          const hsize_t mc_dim_extended[RANK] = 
          {
            offset[0]  +  mc_dim_extend[0],   mc_dim_extend[1],   mc_dim_extend[2]
          };

          // Extend to the new dimension
          ////////hcal_data_set.extend(hcal_dim_extended);
          ////////ecal_data_set.extend(ecal_dim_extended);


          for ( size_t si=0; si<n_subsystems; si++ ) {
             //printf("  %u %s : calling extend with  %u, %u, %u\n", si, subsystem_prefixes[si], subsys_dim_extended[si][0], subsys_dim_extended[si][1], subsys_dim_extended[si][2] ) ;
             calo_data_set[si] -> extend( subsys_dim_extended[si] ) ;
          }

          mc_data_set.extend(mc_dim_extended);

          // Select the hyperslab that only encompasses the
          // difference from extending the data space (i.e. the
          // new event, but offset at the existing event)
          ///////////H5::DataSpace hcal_file_space = hcal_data_set.getSpace();
          ///////////H5::DataSpace ecal_file_space = ecal_data_set.getSpace();

          H5::DataSpace* subsys_file_space[MAX_SUBSYS] ;
          for ( size_t si=0; si<n_subsystems; si++ ) {
             subsys_file_space[si] = new H5::DataSpace( calo_data_set[si] -> getSpace() ) ;
          } // si

          H5::DataSpace mc_file_space = mc_data_set.getSpace();

//        hcal_file_space.selectHyperslab(
//            H5S_SELECT_SET, hcal_dim_extend, offset);

//        ecal_file_space.selectHyperslab(
//            H5S_SELECT_SET, ecal_dim_extend, offset);


          for ( size_t si=0; si<n_subsystems; si++ ) {
             subsys_file_space[si] -> selectHyperslab( H5S_SELECT_SET, calo_dim_extend[si], offset ) ;
          } // si


          mc_file_space.selectHyperslab(
              H5S_SELECT_SET, mc_dim_extend, offset);




          // define memory size to fit the extended hyperslab
          /////////H5::DataSpace hcal_memory_space(RANK, hcal_dim_extend, NULL);
          /////////H5::DataSpace ecal_memory_space(RANK, ecal_dim_extend, NULL);

          H5::DataSpace* subsys_memory_space[MAX_SUBSYS] ;
          for ( size_t si=0; si<n_subsystems; si++ ) {
             subsys_memory_space[si] = new H5::DataSpace( RANK, calo_dim_extend[si], NULL ) ;
          } // si

          H5::DataSpace mc_memory_space(RANK, mc_dim_extend, NULL);





//        // Write the data from memory space to file space
//        hcal_data_set.write(&hcal_data[0], H5::PredType::NATIVE_FLOAT,
//            hcal_memory_space, hcal_file_space);

//        ecal_data_set.write(&ecal_data[0], H5::PredType::NATIVE_FLOAT,
//            ecal_memory_space, ecal_file_space);


          for ( size_t si=0; si<n_subsystems; si++ ) {
             calo_data_set[si] -> write( &(*(subsys_data[si]))[0], H5::PredType::NATIVE_FLOAT,
                  *subsys_memory_space[si], *subsys_file_space[si] ) ;
          } // si

          mc_data_set.write(&mc_data[0], H5::PredType::NATIVE_FLOAT,
              mc_memory_space, mc_file_space);





//        std::fill(hcal_data.begin(),hcal_data.end(),NAN);
//        std::fill(ecal_data.begin(),ecal_data.end(),NAN);

          for ( size_t si=0; si<n_subsystems; si++ ) {
             std::fill( (*(subsys_data[si])).begin(), (*(subsys_data[si])).end(),NAN);
          } // si



          std::fill(mc_data.begin(),mc_data.end(),NAN);


          for ( size_t si=0; si<n_subsystems; si++ ) {
             delete subsys_file_space[si] ;
             delete subsys_memory_space[si] ;
          }

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

    for ( int si=0; si<n_subsystems; si++ ) {
       delete subsysE[si] ;
       delete subsysT[si] ;
       delete subsysX[si] ;
       delete subsysY[si] ;
       delete subsysZ[si] ;
       delete subsys_data[si] ;
    }

  }//argv files
}//write function

//-----------------------------------------------------------------------------------------------------

int main(int argc, char *argv[]){
  if (argc < 3) {
    fprintf(stderr, "%s", "\n\nSyntax is [command] [root_file] [new hdf5 file name]\n\n");
    exit(EXIT_FAILURE);
  }

    char root_file[1000] ;
    sprintf( root_file, "%s", argv[1] ) ;

    printf("\n\n Opening input root file %s\n\n", root_file ) ;

    TFile *tf_file = TFile::Open(root_file);
    if ( tf_file == 0x0 ) { printf("\n\n *** Bad input file:  %s\n\n", root_file ) ; exit(-1) ; }
    if ( ! tf_file -> IsOpen() ) { printf("\n\n *** Bad input file:  %s\n\n", root_file ) ; exit(-1) ; }


    TTree* tt_events = (TTree*) tf_file -> Get( "events" ) ;
    if ( tt_events == 0x0 ) { printf("\n\n *** can't find events TTree in %s\n\n", root_file ) ; exit(-1) ; }

    TObjArray* lob = tt_events -> GetListOfBranches() ;
    printf("\n\n Branches in input file:\n") ;
    {
       int si=0; 
       bool found_mcparticles = false ;
       for ( int i=0; i< lob->GetEntries(); i++ ) {
          printf( "  %3d : %s\n", i, lob->At(i) -> GetName() ) ;
          char ssname[1000] ;
          sprintf( ssname, "HcalEndcapPHitsReco" ) ;
          if ( strcmp( lob->At(i) -> GetName(), ssname ) == 0 ) {
             printf("             - found %s.\n", ssname ) ;
             sprintf( subsystem_prefixes[si], "%s", ssname ) ;
             sprintf( subsystem_short_names[si], "hcal" )  ;
             si++ ;
          }
          sprintf( ssname, "HcalEndcapPInsertHitsReco" ) ;
          if ( strcmp( lob->At(i) -> GetName(), ssname ) == 0 ) {
             printf("             - found %s.\n", ssname ) ;
             sprintf( subsystem_prefixes[si], "%s", ssname ) ;
             sprintf( subsystem_short_names[si], "hcali" )  ;
             si++ ;
          }
          sprintf( ssname, "EcalEndcapPHitsReco" ) ;
          if ( strcmp( lob->At(i) -> GetName(), ssname ) == 0 ) {
             printf("             - found %s.\n", ssname ) ;
             sprintf( subsystem_prefixes[si], "%s", ssname ) ;
             sprintf( subsystem_short_names[si], "ecal" )  ;
             si++ ;
          }
          sprintf( ssname, "EcalEndcapPInsertHitsReco" ) ;
          if ( strcmp( lob->At(i) -> GetName(), ssname ) == 0 ) {
             printf("             - found %s.\n", ssname ) ;
             sprintf( subsystem_prefixes[si], "%s", ssname ) ;
             sprintf( subsystem_short_names[si], "ecali" )  ;
             si++ ;
          }
          if ( strcmp( lob->At(i) -> GetName(), "MCParticles" ) == 0 ) {
             found_mcparticles = true ;
             printf("             - found MCParticles.\n", ssname ) ;
          }
       }
       n_subsystems = si ;
       if ( !found_mcparticles ) { printf("\n\n *** Didn't find MCParticles in events TTree.\n\n") ; exit(-1) ; }
    }
    printf("\n\n found %d subsystems.\n\n", n_subsystems ) ;



  static const size_t cal_row_size = 4; //Number of calorimeter hit variables
  static const size_t mc_row_size = 10; //Number of MC truth variables

  size_t eventsN_max = 0;
  size_t calo_NHits_max = 0;
  size_t mcNParticles_max = 0;
  size_t block_size = 100; //affects chunk size, 

  //Constants for Layering HCal
  const double z_offset = 3800.; //[mm]. Fixes some hardcoded setting in ATHENA detector
  const double z_max = 1200.; // actual length of hcal in z [mm]

  find_max_dims(argv + 1, argv + argc - 1, eventsN_max, calo_NHits_max, mcNParticles_max);
  /* eventsN_max = 10000; calo_NHits_max = 1318*2; mcNParticles_max = 15; */ 
  //saved for rec_piplus_Energy_0-100GeV.root





  // Access mode H5F_ACC_TRUNC truncates any existing file, while
  // not throwing an exception (unlike H5F_ACC_RDWR)
  std::string file_str = argv[argc-1];
  H5::H5File file( file_str.c_str(), H5F_ACC_TRUNC );
  // The tensor dimension for each new chunk of events

  //The chunking of data can be edited for performance
  ////////hsize_t hcal_dim_extend[RANK] = {block_size,  cal_row_size, calo_NHits_max};
  ////////hsize_t ecal_dim_extend[RANK] = {block_size, cal_row_size, calo_NHits_max};

  hsize_t calo_dim_extend[n_subsystems][RANK] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
     calo_dim_extend[si][0] = block_size ;
     calo_dim_extend[si][1] = cal_row_size ;
     calo_dim_extend[si][2] = calo_NHits_max ;
  } // si


  hsize_t mc_dim_extend[RANK] = {block_size, mc_row_size, mcNParticles_max};

  //Check the hyperslab/extension dimensions are correct
  fprintf(stderr,"\n%s: %d: HDF5 Chunk Size = %u\n",__func__,__LINE__,block_size);

  for ( size_t si=0; si<n_subsystems; si++ ) {
     fprintf(stderr, "%s: %d: %s : calo_dim_extend = %u %u %u\n", 
      __func__, __LINE__, subsystem_prefixes[si], calo_dim_extend[si][0],calo_dim_extend[si][1],calo_dim_extend[si][2]);
  }
  fprintf(stderr, "%s: %d:  mc_dim_extend  = %u %u  %u\n", 
      __func__, __LINE__, mc_dim_extend[0],mc_dim_extend[1],mc_dim_extend[2]);


  // The maximum tensor dimension, for unlimited number of events
  // a.k.a. the overall dimensions of the dataset
  /////////////hsize_t hcal_dim_max[RANK] = {H5S_UNLIMITED, cal_row_size, calo_NHits_max};
  /////////////hsize_t ecal_dim_max[RANK] = {H5S_UNLIMITED, cal_row_size, calo_NHits_max};

  hsize_t calo_dim_max[n_subsystems][RANK] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
     calo_dim_max[si][0] = H5S_UNLIMITED ;
     calo_dim_max[si][1] = cal_row_size ;
     calo_dim_max[si][2] = calo_NHits_max ;
  } // si




  hsize_t mc_dim_max[RANK] =  {H5S_UNLIMITED, mc_row_size, mcNParticles_max};

  // The extensible HDF5 data space
  //////////H5::DataSpace hcal_data_space(RANK, hcal_dim_extend, hcal_dim_max);
  //////////H5::DataSpace ecal_data_space(RANK, ecal_dim_extend, ecal_dim_max);

  H5::DataSpace* calo_data_space[MAX_SUBSYS] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
     printf("  Making dataspace for subsystem %d %s\n", si, subsystem_prefixes[si] ) ;
     calo_data_space[si] = new H5::DataSpace( RANK, calo_dim_extend[si], calo_dim_max[si] ) ;
  } // si




  printf("  Making dataspace for MCParticles\n" ) ;
  H5::DataSpace mc_data_space(RANK, mc_dim_extend, mc_dim_max);

  // To enable zlib compression (there will be many NANs) and
  // efficient chunking (splitting of the tensor into contingous
  // hyperslabs), a HDF5 property list is needed
  ///////////H5::DSetCreatPropList hcal_property = H5::DSetCreatPropList();
  ///////////H5::DSetCreatPropList ecal_property = H5::DSetCreatPropList();

  H5::DSetCreatPropList* calo_property[MAX_SUBSYS] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
     calo_property[si] = new H5::DSetCreatPropList();
  } // si

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
      ////////hcal_property.setDeflate(1);
      ////////ecal_property.setDeflate(1);

      printf(" calling setDeflate\n") ;
      for ( size_t si=0; si<n_subsystems; si++ ) {
         calo_property[si] -> setDeflate(1) ;
      }
      mc_property.setDeflate(1);
    }
  }
#endif // HDF5_USE_DEFLATE

  // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
  // the CPU L2 cache size

  /////////hsize_t hcal_dim_chunk[RANK] = {
  /////////  hcal_dim_extend[0],
  /////////  hcal_dim_extend[1],
  /////////  hcal_dim_extend[2]
  /////////};

  /////////hsize_t ecal_dim_chunk[RANK] = {
  /////////  ecal_dim_extend[0],
  /////////  ecal_dim_extend[1],
  /////////  ecal_dim_extend[2]
  /////////};

  hsize_t calo_dim_chunk[MAX_SUBSYS][RANK] ;
  for ( size_t si; si<n_subsystems; si++ ) {
     calo_dim_chunk[si][0] = calo_dim_extend[si][0] ;
     calo_dim_chunk[si][1] = calo_dim_extend[si][1] ;
     calo_dim_chunk[si][2] = calo_dim_extend[si][2] ;
  }

  hsize_t mc_dim_chunk[RANK] = {
    mc_dim_extend[0],
    mc_dim_extend[1],
    mc_dim_extend[2]
  };

  ///////hcal_property.setChunk(RANK, hcal_dim_chunk);
  ///////ecal_property.setChunk(RANK, ecal_dim_chunk);
  for ( size_t si; si<n_subsystems; si++ ) {
     printf("  Calling setChunk for subsystem %u  %s :  %u, %u, %u\n", si, subsystem_prefixes[si], calo_dim_chunk[si][0], calo_dim_chunk[si][1], calo_dim_chunk[si][2] ) ;
     calo_property[si] -> setChunk( RANK, calo_dim_chunk[si] ) ;
  } // si

  mc_property.setChunk(RANK, mc_dim_chunk);


  // Create the data set, which will have space for the first event (chunk)
//H5::DataSet hcal_data_set =
//  file.createDataSet("hcal", H5::PredType::NATIVE_FLOAT,
//      hcal_data_space, hcal_property);

//H5::DataSet ecal_data_set =
//  file.createDataSet("ecal", H5::PredType::NATIVE_FLOAT,
//      ecal_data_space, ecal_property);


  H5::DataSet* calo_data_set[MAX_SUBSYS] ;
  for ( size_t si; si<n_subsystems; si++ ) {
     calo_data_set[si] = new H5::DataSet( file.createDataSet( subsystem_short_names[si], H5::PredType::NATIVE_FLOAT, *calo_data_space[si], *calo_property[si] ) ) ;
  } // si


  H5::DataSet mc_data_set =
    file.createDataSet("mc", H5::PredType::NATIVE_FLOAT,
        mc_data_space, mc_property);

  hsize_t offset[RANK] = {0, 0, 0};

//write_data(hcal_data_set, ecal_data_set, mc_data_set,
//    cal_row_size, mc_row_size,
//    offset, hcal_dim_extend, ecal_dim_extend, mc_dim_extend,
//    eventsN_max, calo_NHits_max, mcNParticles_max, block_size,argv + 1, argv + argc - 1);

  write_data(calo_data_set, mc_data_set,
      cal_row_size, mc_row_size,
      offset, calo_dim_extend, mc_dim_extend,
      eventsN_max, calo_NHits_max, mcNParticles_max, block_size,argv + 1, argv + argc - 1);




  fprintf(stderr,"\n\n%s: %d: [Complete] \n\n",__func__,__LINE__);

  file.close();

  for ( size_t si; si<n_subsystems; si++ ) {
     delete calo_data_set[si] ;
     delete calo_data_space[si] ;
     delete calo_property[si] ;
  } // si

  return EXIT_SUCCESS;
}




