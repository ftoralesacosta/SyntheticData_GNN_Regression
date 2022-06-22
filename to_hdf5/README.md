## This directory contains code for converting reconstruction root files generated using `eic/reconstruction_benchmarks`.

### Convert the root file to hdf5:

> ./root_to_hdf5 [reconstruction_root_file] [new event_hdf5_file]

The data of the resulting `event hdf5 file` are containde in `mc`, `hcal` and `ecal` dataspaces. The data in those dataspaces is as follows:

|     |                    |  Generator Data:  |      |       |       |      |       |     |          |
|:---:|:------------------:|:-----------------:|:----:|:-----:|:-----:|:----:|:-----:|:---:|:--------:|
| PDG | mcSimulatorStatus | mcGeneratorStatus | $P_X$ | $P_Y$ | $P_Z$ | Mass | $P_T$ | $P$ | $\theta$ |
| 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |

|        | ECal:  |        |        |     |        | HCal:  |        |        |
|--------|--------|--------|--------|-----|--------|--------|--------|--------|
| Cell E | Cell X | Cell Y | Cell Z |     | Cell E | Cell X | Cell Y | Cell Z |
| 0      | 1      | 2      | 3      |     | 0      | 1      | 2      | 3      |

---

### Process the hdf5 file for use in training a Neural Network from [this repo](https://github.com/eiccodesign/regressiononly).
##### `events`->`calorimeter images`

> ./H5_GetImages [event_hdf5_file] [new image_hdf5_file]


The data in the resulting `image hdf5 file` are in `calo_images` and `truth` dataspaces:
|        |        |        |  Calo Image  |           |           |
|:------:|:------:|:------:|:------------:|:---------:|:---------:|
| Cell E | Cell X | Cell Y | Cell 'Depth' | Layer 1 Z | Layer 2 Z |
|    0   |    1   |    2   |       3      |     4     |     5     |

| Truth   |          |
|:-------:|:--------:|
|Truth $E$| $\theta$ |
|    0    |     1    |

__Tip:__ The H5_hitQA notebook is run on the event_hdf5_files, while the H5_TowardsCodesign notebook is run on the image_hdf5_files.
