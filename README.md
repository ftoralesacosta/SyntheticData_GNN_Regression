This Branch will combine several of the repositories listed below, as well as HDF5 files for the final data format.

> ./get_eic-container.sh
This will download the container, enter the container, and create some required directories

> source get_frameworks.sh
This downloads the specific commits used to generate data. It builds them, and sources the setup_env.sh script to set the proper paths
Make sure you're still in the container when running this.

> bash benchmarks/clustering/full_cal_clusters.sh -p "electron" -n 100 --pmin 0.99 --pmax 1.01 -t test1
Also make sure to still be in the container when running this.
This last command will use single-particle gun to generate 100 electrons with a momentum in range 0.99-1.01 GeV. 
Then, digitization and reconstruction will be run. 
Output will contain sim_TAGNAME.root and rec_TAGNAME.root files, which correspond to the G4-level and reconstructed level respectively (digitization, clustering). 
The files will contain a ROOT TTree with the generated particles, hits at the G4 level, hits after digitization, and clusters. 
We will mainly be working with hits either before or after digitization.

For HDF5: While still inside the container, navigate to `to_hdf5`, and execute the grab hdf5 script. May need more testing for portability as of 4/28/22
TODO: Test portabiliby of HDF5
