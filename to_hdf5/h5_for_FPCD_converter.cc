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

#define particle_gun true
#define RANK 3 //Event No., branch, variable index
#define cluster_RANK 2 //Event No., branch, variable index

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

//Obtain No. cell hits and MC particles for setting HDF5 dimensions
void find_max_dims(
    char *argv_first[],
    char *argv_last[],
    size_t &eventsN_max,
    size_t &cell_NHits_max,
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
               cell_NHits_max = std::max(cell_NHits_max, subsysNHits);
          }

          if (i%1000 == 0) {
            fprintf(stderr, "%s: %d: event = %i / %u, max cell Hits= %u, max MC Particles = %u\n",
                __func__, __LINE__, i, eventsN_max, cell_NHits_max, mcNParticles_max);
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

    printf("\n\n  find_max_dims :   eventsN_max = %u,  cell_NHits_max = %u,  mcNParticles_max = %u\n\n",
       eventsN_max, cell_NHits_max, mcNParticles_max ) ;

  }//argc

}//find_max_dims

//-----------------------------------------------------------------------------------------------------

void write_data(
    H5::DataSet* cell_data_set[],
    H5::DataSet &cluster_data_set,
    size_t cell_row_size,
    size_t cluster_row_size,
    hsize_t *offset,
    hsize_t *cluster_offset,
    const hsize_t cell_dim_extend[][RANK],
    const hsize_t *cluster_dim_extend,
    const UInt_t eventsN_max,
    const UInt_t cell_NHits_max,
    const UInt_t mcNParticles_max,
    const UInt_t block_size,
    char *argv_first[], char *argv_last[]) 
{

  printf("\n\n ========== Beginning of write_data\n\n") ;

  for (char **p = argv_first; p != argv_last; p++) {

    TFile *file = TFile::Open(*p);
    TTreeReader events("events", file);

    char mc_str[]  =  "MCParticles";

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


    std::vector<float>* subsys_data[MAX_SUBSYS] ;
    for ( int si=0; si<n_subsystems; si++ ) {
       subsys_data[si] = new std::vector<float>( block_size * cell_row_size * cell_NHits_max, 0 );
    } // si


    std::vector<float> cluster_data(block_size *  cluster_row_size, 0);

    //Check the max dims are passed correctly

    int i = 0;
    while (events.Next()) {

      if (i >= eventsN_max) break;
      int iblock = i % block_size;
      //writing to file is done every [block_size] number of events
      //this variable keeps track of the current increment within a block,
      //as opposed to [i] which is looping through all events. Note 'i' is 
      //incremented only after an event passes selection (hit requirements)

      size_t subsys_fill[MAX_SUBSYS] ;
      float_t subsys_sumE[MAX_SUBSYS] ;
      for ( int si=0; si<n_subsystems; si++ ) 
      { 
        subsys_fill[si] = 0 ; 
        subsys_sumE[si] = 0;
      }
      //counts cell hits that pass cuts

      size_t subsysNHits[MAX_SUBSYS] ;
      for ( int si=0; si<n_subsystems; si++ ) { subsysNHits[si] = subsysE[si] -> GetSize() ; }
      //Gets Max N hits, helpfull for looping and setting array sizes

      for ( int si=0; si<n_subsystems; si++ ) {
        for (size_t h_hit = 0; h_hit < subsysNHits[si]; h_hit++) {
            if ( (*(subsysE[si]))[h_hit] > 1e10) continue ; //Spikey Cells. Should be fixed in Juggler Commit post Feb2022
            if ( (*(subsysE[si]))[h_hit] <= 0.00006 ) continue ; //MIPS and Empty Cells
            if ( (*(subsysT[si]))[h_hit] > 200 ) continue ;

            size_t E_index = iblock*cell_NHits_max*cell_row_size + (cell_row_size*subsys_fill[si]) + 0;
            size_t X_index = iblock*cell_NHits_max*cell_row_size + (cell_row_size*subsys_fill[si]) + 1;
            size_t Y_index = iblock*cell_NHits_max*cell_row_size + (cell_row_size*subsys_fill[si]) + 2;
            size_t Z_index = iblock*cell_NHits_max*cell_row_size + (cell_row_size*subsys_fill[si]) + 3;
            size_t Mask_index=iblock*cell_NHits_max*cell_row_size+ (cell_row_size*subsys_fill[si]) + 4;
            //Index for flattened 3D vector

            (*(subsys_data[si]))[E_index] = (*(subsysE[si]))[h_hit] ; //GeV
            (*(subsys_data[si]))[X_index] = (*(subsysX[si]))[h_hit] ;
            (*(subsys_data[si]))[Y_index] = (*(subsysY[si]))[h_hit] ;
            (*(subsys_data[si]))[Z_index] = (*(subsysZ[si]))[h_hit] ;
            (*(subsys_data[si]))[Mask_index] = 1;//See JetNet - arXiv:2106.11535 "mask"

            subsys_fill[si] ++ ;
            subsys_sumE[si] += (*(subsys_data[si]))[E_index];

         } // h_hit

         /* fprintf(stderr, "%u %s: n_fill = %u \n",__LINE__,__func__,subsys_fill[si]); */
         /* fprintf(stderr, "%u %s: sumE = %1.4f \n",__LINE__,__func__,subsys_sumE[si]); */
         if ( subsysNHits[si] == 0 ) continue ;
         if ( subsys_fill[si] == 0 ) continue ;
         if ( subsys_sumE[si]  < 0.1   ) continue ;

      }//si

      size_t mc_fill = 0;
      size_t mcNParticles = mcMass.GetSize();
      for (size_t particle = 0; particle < mcNParticles; particle++) 
      {
        if (mcGeneratorStatus[particle] != 1) continue;

        //Calculated Values
        float mcPT = hypot(mcPX[particle], mcPY[particle]);
        float mcP = hypot(mcPX[particle], mcPY[particle], mcPZ[particle]);
        float mcTheta = acos(mcPZ[particle]/mcP)*180./M_PI;

        cluster_data[iblock*cluster_row_size + 0] = mcP;
        cluster_data[iblock*cluster_row_size + 1] = mcTheta;
        /* mc_data[(iblock*mc_row_size + 0)*mcNParticles_max + mc_fill] = mcPDG[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 1)*mcNParticles_max + mc_fill] = mcSimulatorStatus[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 2)*mcNParticles_max + mc_fill] = mcGeneratorStatus[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 3)*mcNParticles_max + mc_fill] = mcPX[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 4)*mcNParticles_max + mc_fill] = mcPY[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 5)*mcNParticles_max + mc_fill] = mcPZ[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 6)*mcNParticles_max + mc_fill] = mcMass[particle]; */ 
        /* mc_data[(iblock*mc_row_size + 7)*mcNParticles_max + mc_fill] = mcPT; */ 
        if (particle_gun)
          break;//break after first particle
      }

      //FIXME: Fill with MCParticles Momentum
      cluster_data[iblock*cluster_row_size + 2] = subsys_sumE[0]; 
      cluster_data[iblock*cluster_row_size + 3] = subsys_fill[0]; //FIXME: ONLY READS FIRST SUBSYSTEM, make sure HCAL
      //FIXME: Put outside fo se look, chang to  Nsub*row size

      bool print_cal = false;
      bool print_mc = false;


      if (iblock == (block_size-1)) 
      {
        //writes 1 block (100 events) at a time. Faster/less memory

        if (print_cal){
          for ( int si=0; si<n_subsystems; si++ ) {
            for (size_t h_hit = 0; h_hit < subsysNHits[si]; h_hit++) {
              float E = (*(subsys_data[si]))[(iblock*cell_row_size + 0)*cell_NHits_max + h_hit];
              float Z = (*(subsys_data[si]))[(iblock*cell_row_size + 1)*cell_NHits_max + h_hit];
              if (std::isnan(E)) break;
              if (h_hit%10 == 0)
              {
                printf("   %s Hit # %5d / %5d, E = %10.6f\n", subsystem_prefixes[si], h_hit, subsysNHits[si], E ) ;
                printf("   %s Hit # %5d / %5d, Z = %10.6f\n", subsystem_prefixes[si], h_hit, subsysNHits[si], Z ) ;
              }
            } // h_hit
          } // si
        } // print_cal?

        if (print_mc) {
          fprintf(stderr, "\n%u %s: Number of Cells = %llu\n",
              __LINE__,__func__,cluster_data[iblock*cluster_row_size + 0]);
          fprintf(stderr, "%u %s: Cluster Sum = %1.2f\n",
              __LINE__,__func__,cluster_data[iblock*cluster_row_size + 1]);
        }//print

        if (offset[0] == 0) {
          // Writing the first event. Data spaces is already created with space for one event
          // see file.createDataSet()
          for ( size_t si=0; si<n_subsystems; si++ ) {
            cell_data_set[si] -> write( &(*(subsys_data[si]))[0], H5::PredType::NATIVE_FLOAT );
          } // si
          cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT);

          //Make sure to clear previous arrays
          for ( size_t si=0; si<n_subsystems; si++ ) {
            std::fill( (*(subsys_data[si])).begin(), (*(subsys_data[si])).end(),0);
          }
          std::fill(cluster_data.begin(),cluster_data.end(),0);
        }

        else { //offset[0] =/= 0

          hsize_t subsys_dim_extended[MAX_SUBSYS][RANK] ;

          for ( size_t si=0; si<n_subsystems; si++ ) {
            subsys_dim_extended[si][0] = offset[0] + cell_dim_extend[si][0] ;
            subsys_dim_extended[si][1] = cell_dim_extend[si][1] ;
            subsys_dim_extended[si][2] = cell_dim_extend[si][2] ;
          } // si


          const hsize_t cluster_dim_extended[cluster_RANK] = 
          { cluster_offset[0]  +  cluster_dim_extend[0],   cluster_dim_extend[1] };

          // Extend to the new dimension
          for ( size_t si=0; si<n_subsystems; si++ ) 
          {
            cell_data_set[si] -> extend( subsys_dim_extended[si] ) ;
          }

          cluster_data_set.extend(cluster_dim_extended);

          // Select the hyperslab that only encompasses the
          // difference from extending the data space (i.e. the
          // new event, but offset at the existing event)

          H5::DataSpace* subsys_file_space[MAX_SUBSYS] ;
          for ( size_t si=0; si<n_subsystems; si++ ) {
            subsys_file_space[si] = new H5::DataSpace( cell_data_set[si] -> getSpace() ) ;
          } // si

          for ( size_t si=0; si<n_subsystems; si++ ) {
            subsys_file_space[si] -> selectHyperslab( H5S_SELECT_SET, cell_dim_extend[si], offset ) ;
          } // si

          H5::DataSpace cluster_file_space = cluster_data_set.getSpace();

          cluster_file_space.selectHyperslab(
              H5S_SELECT_SET, cluster_dim_extend, cluster_offset);


          // define memory size to fit the extended hyperslab
          H5::DataSpace* subsys_memory_space[MAX_SUBSYS] ;
          for ( size_t si=0; si<n_subsystems; si++ ) {
            subsys_memory_space[si] = new H5::DataSpace( RANK, cell_dim_extend[si], NULL ) ;
          } // si
          for ( size_t si=0; si<n_subsystems; si++ ) {
            cell_data_set[si] -> write( &(*(subsys_data[si]))[0], H5::PredType::NATIVE_FLOAT,
                *subsys_memory_space[si], *subsys_file_space[si] ) ;
          } // si

          H5::DataSpace cluster_memory_space(cluster_RANK, cluster_dim_extend, NULL);
          cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT,
              cluster_memory_space, cluster_file_space);


          // Remember to clear arrays before beginning next block!
          for ( size_t si=0; si<n_subsystems; si++ ) {
            std::fill( (*(subsys_data[si])).begin(), (*(subsys_data[si])).end(),0);
          } // si

          std::fill(cluster_data.begin(),cluster_data.end(),0);

          //Delete H5 file and memory spaces
          for ( size_t si=0; si<n_subsystems; si++ ) {
            delete subsys_file_space[si] ;
            delete subsys_memory_space[si] ;
          }

        }//Offset[0] =/= 0? 

        offset[0]+=block_size;
        cluster_offset[0]+=block_size;

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
        sprintf( subsystem_short_names[si], "hcal_cells" )  ;
        si++ ;
      }
      sprintf( ssname, "HcalEndcapPInsertHitsReco" ) ;
      if ( strcmp( lob->At(i) -> GetName(), ssname ) == 0 ) {
        printf("             - found %s.\n", ssname ) ;
        printf("             - NOT USING, SKIPPING %s.\n", ssname ) ;
        continue;
        /* sprintf( subsystem_prefixes[si], "%s", ssname ) ; */
        /* sprintf( subsystem_short_names[si], "hcali" )  ; */
        /* si++ ; */
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


  static const size_t cell_row_size = 5; //Number of cellrimeter hit variables
  static const size_t cluster_row_size = 4; //Cluster Variables: Energy, N_Cells

  size_t eventsN_max = 0;
  size_t cell_NHits_max = 0;
  size_t mcNParticles_max = 0;
  size_t block_size = 100; //affects chunk size, 

  //Constants for Layering HCal
  const double z_offset = 3800.; //[mm]. Fixes some hardcoded setting in ATHENA detector
  const double z_max = 1200.; // actual length of hcal in z [mm]

  find_max_dims(argv + 1, argv + argc - 1, eventsN_max, cell_NHits_max, mcNParticles_max);
  /* eventsN_max = 200000; cell_NHits_max = 1861; mcNParticles_max = 30; */


  // Access mode H5F_ACC_TRUNC truncates any existing file, while
  // not throwing an exception (unlike H5F_ACC_RDWR)
  std::string file_str = argv[argc-1];
  H5::H5File file( file_str.c_str(), H5F_ACC_TRUNC );

  //The tensor dimension for each new chunk of events
  hsize_t cell_dim_extend[n_subsystems][RANK] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
    cell_dim_extend[si][0] = block_size ;
    cell_dim_extend[si][1] = cell_NHits_max ;
    cell_dim_extend[si][2] = cell_row_size ;
  } // si

  hsize_t cluster_dim_extend[cluster_RANK] = {block_size, cluster_row_size};

  //Check the hyperslab/extension dimensions are correct
  fprintf(stderr,"\n%s: %d: HDF5 Chunk Size = %u\n",__func__,__LINE__,block_size);

  for ( size_t si=0; si<n_subsystems; si++ ) {
    fprintf(stderr, "%s: %d: %s : cell_dim_extend = %u %u %u\n", 
        __func__, __LINE__, subsystem_prefixes[si], cell_dim_extend[si][0],cell_dim_extend[si][1],cell_dim_extend[si][2]);
  }
  fprintf(stderr, "%s: %d:  cluster_dim_extend  = %u %u  %u\n", 
      __func__, __LINE__, cluster_dim_extend[0],cluster_dim_extend[1]);


  //Create cell Subsystem Dataspaces + PropLists
  hsize_t cell_dim_max[n_subsystems][RANK] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
    cell_dim_max[si][0] = H5S_UNLIMITED ; //for unlimited NEvents
    cell_dim_max[si][1] = cell_NHits_max ;
    cell_dim_max[si][2] = cell_row_size ; //N cell variables per hit
  } // si

  H5::DataSpace* cell_data_space[MAX_SUBSYS] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
    printf("  Making dataspace for subsystem %d %s\n", si, subsystem_prefixes[si] ) ;
    cell_data_space[si] = new H5::DataSpace( RANK, cell_dim_extend[si], cell_dim_max[si] ) ;
  } // si

  H5::DSetCreatPropList* cell_property[MAX_SUBSYS] ;
  for ( size_t si=0; si<n_subsystems; si++ ) {
    cell_property[si] = new H5::DSetCreatPropList();
  } // si


  //Create MC Gen Dataspaces + PropLists
  hsize_t cluster_dim_max[cluster_RANK] =  {H5S_UNLIMITED, cluster_row_size};
  printf("  Making dataspace for MCParticles\n" ) ;
  H5::DataSpace cluster_data_space(cluster_RANK, cluster_dim_extend, cluster_dim_max);
  H5::DSetCreatPropList cluster_property = H5::DSetCreatPropList();

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
      printf(" calling setDeflate\n") ;
      for ( size_t si=0; si<n_subsystems; si++ ) {
        cell_property[si] -> setDeflate(1) ;
      }
      cluster_property.setDeflate(1);
    }
  }
#endif // HDF5_USE_DEFLATE

  hsize_t cell_dim_chunk[MAX_SUBSYS][RANK] ;
  for ( size_t si; si<n_subsystems; si++ ) {
    cell_dim_chunk[si][0] = cell_dim_extend[si][0] ;
    cell_dim_chunk[si][1] = cell_dim_extend[si][1] ;
    cell_dim_chunk[si][2] = cell_dim_extend[si][2] ;
  }

  hsize_t cluster_dim_chunk[cluster_RANK] = {
    cluster_dim_extend[0],
    cluster_dim_extend[1]
  };


  for ( size_t si; si<n_subsystems; si++ ) {
    printf("  Calling setChunk for subsystem %u  %s :  %u, %u, %u\n", si, subsystem_prefixes[si], cell_dim_chunk[si][0], cell_dim_chunk[si][1], cell_dim_chunk[si][2] ) ;
    cell_property[si] -> setChunk( RANK, cell_dim_chunk[si] ) ;
  } // si

  cluster_property.setChunk(cluster_RANK, cluster_dim_chunk);


  H5::DataSet* cell_data_set[MAX_SUBSYS] ;
  for ( size_t si; si<n_subsystems; si++ ) {
    cell_data_set[si] = new H5::DataSet( file.createDataSet( subsystem_short_names[si], H5::PredType::NATIVE_FLOAT, *cell_data_space[si], *cell_property[si] ) ) ;
  } // si


  H5::DataSet cluster_data_set =
    /* file.createDataSet("mc", H5::PredType::NATIVE_FLOAT, */
    file.createDataSet("cluster", H5::PredType::NATIVE_FLOAT,
        cluster_data_space, cluster_property);

  hsize_t offset[RANK] = {0, 0, 0};
  hsize_t cluster_offset[cluster_RANK] = {0, 0};

  write_data(cell_data_set, cluster_data_set,
      cell_row_size, cluster_row_size, offset,  cluster_offset, 
      cell_dim_extend, cluster_dim_extend,
      eventsN_max, cell_NHits_max, mcNParticles_max, 
      block_size,argv + 1, argv + argc - 1);

  fprintf(stderr,"\n\n%s: %d: [Complete] \n\n",__func__,__LINE__);

  file.close();

  for ( size_t si; si<n_subsystems; si++ ) {
    delete cell_data_set[si] ;
    delete cell_data_space[si] ;
    delete cell_property[si] ;
  } // si

  return EXIT_SUCCESS;
}

