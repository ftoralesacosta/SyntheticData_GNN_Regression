### This directory contains code for converting reconstruction root files generated using `eic/reconstruction_benchmarks`.

#### Convert the root file to hdf5:

> ./root_to_hdf5 [reconstruction_root_file] [new event_hdf5_file]

#### Process the hdf5 file for use in training a Neural Network from [this repo](https://github.com/eiccodesign/regressiononly).
##### `events`->`calorimeter images`

> ./H5_GetImages [event_hdf5_file] [new image_hdf5_file]

__Tip:__ The H5_hitQA notebook is run on the event_hdf5_files, while the H5_TowardsCodesign notebook is run on the image_hdf5_files.
