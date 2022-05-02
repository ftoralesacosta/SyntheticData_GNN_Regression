CXX =	g++
CC = $(CXX)
CXXFLAGS = -g -O -s -frounding-math
CXXFLAGS += -std=c++17

CXXFLAGS +=	$(shell root-config --cflags)
LDFLAGS =	$(shell root-config --ldflags)
LDLIBS =	$(shell root-config --libs)

HDF5_ROOT=$ENV{HDF5_DIR}
$(info $$HDF5_DIR is [${HDF5_DIR}])
ifneq (,$(HDF5_DIR))
    CXXFLAGS +=	-I$(HDF5_DIR)/include
    LDFLAGS +=	-L$(HDF5_DIR)/lib64 -L$(HDF5_DIR)/lib -L$(HDF5_DIR)
    LDLIBS +=	-lhdf5_cpp -lhdf5
    #LDLIBS +=	-lhdf5
endif

# include_directories(${HDF5_INCLUDE_DIRS})

CXXFLAGS +=	-I../.. -I.

TARGET = root_to_hdf5	


