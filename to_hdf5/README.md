## This directory contains code for converting reconstruction root files generated using `eic/reconstruction_benchmarks`.

### Compile the code:

The makefile works on both __root_to_hdf5.cc__ and __H5_GetImages.cc__, one simply has to run __make__ and specify the target

> make root_to_hdf5

> make H5_GetImages

### Convert the root file to hdf5:

> ./root_to_hdf5 [reconstruction_root_file] [new event_hdf5_file]

The data of the resulting `event hdf5 file` are contained in __mc__, __hcal__ and __ecal__ dataspaces. The data in those dataspaces is structured like this:

|     |                    |  Generator Data:  |      |       |       |      |       |     |          |
|:---:|:------------------:|:-----------------:|:----:|:-----:|:-----:|:----:|:-----:|:---:|:--------:|
| PDG | mcSimulatorStatus  | mcGeneratorStatus | $P_X$| $P_Y$ | $P_Z$ | Mass | $P_T$ | $P$ | $\theta$ |
|  0  |         1         |          2         |   3  |   4   |   5   |   6  |   7   |  8  |    9     |

|        | ECal:  |        |        |     |        | HCal:  |        |        |
|--------|--------|--------|--------|-----|--------|--------|--------|--------|
| Cell E | Cell X | Cell Y | Cell Z |     | Cell E | Cell X | Cell Y | Cell Z |
| 0      | 1      | 2      | 3      |     | 0      | 1      | 2      | 3      |

---

### Process the hdf5 file for use in training a Neural Network from [this repo](https://github.com/eiccodesign/regressiononly).
#### `events` $\rightarrow$`calorimeter images`

> ./H5_GetImages [event_hdf5_file] [new image_hdf5_file]


The data in the resulting `image hdf5 file` are in __calo_images__ and __truth__ dataspaces. The data is `normalized` according to:
$$ (X_{h,i} - \mu_i)/\sigma_i $$

Where X is the set  input data, index *h* is the cell hit index, and **i** the index corresponding to cell variable. $\ \mu_i\ \mathrm{and}\ \sigma_i$ are the mean and standard deviation of cell variable *i*. The cell variables are structured acording to this table:

|        |        |        |  Calo Image  |           |           |
|:------:|:------:|:------:|:------------:|:---------:|:---------:|
| Cell E | Cell X | Cell Y | Cell 'Depth' | Layer 1 Z | Layer 2 Z |
|    0   |    1   |    2   |       3      |     4     |     5     |

And the **truth** variables are:

| Truth   |          |
|:-------:|:--------:|
|Truth $E$| Truth $\theta$ |
|    0    |     1    |

__Tip:__ The H5_hitQA notebook is run on the event_hdf5_files, while the H5_TowardsCodesign notebook is run on the image_hdf5_files.
