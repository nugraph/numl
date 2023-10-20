# numl larsoft package

This repository contains a LArSoft package for generating numl-format event HDF5 files from artroot files. These event HDF5 files contain low-level information (simulated particles, energy deposits, detector hits etc), which can then be used for efficient downstream production of machine learning inputs using the [pynuml](https://github.com/vhewes/pynuml) package.

## Installation

The `develop` branch of this repository is kept up-to-date with the DUNE flavour of LArSoft, with tagged versions produced regularly. If you have a `dunesw` development area set up via MRB, you can install this package simply by running

```
cd $MRB_SOURCE
mrb g -t <larsoft version tag> https://github.com/vhewes/numl
mrbsetenv
mrb i -j4
```

## Event HDF5 Generation

Once the repository is built, you should be able to run

```
lar -c hdf5maker_dune.fcl <artroot file>
```

to generate an event HDF5 file for downstream processing.

## Event HDF5 Merging

The [pynuml](https://github.com/vhewes/pynuml) package is designed to operate on a single HDF5 file that contains the full dataset. However, since high-statistics MC generation is typically performed at scale via grid computing, the output will typically be a larger number of small event HDF5 files. These files can be merged into a single file using the [ph5concat](https://github.com/NU-CUCIS/ph5concat) package. This package can be installed manually, but for those working on Linux architectures, it is recommended to install via Anaconda as follows:

```
conda install -c conda-forge ph5concat
```

This concatenation utility leverages MPI to merge files efficiently in parallel on High-Performance Computing (HPC) nodes. The user can run

```
# merge files
ph5_concat -n <N> -i files.txt -o out.evt.h5

# add sequencing metadata to files
add_key -k event_table/event_id -c -f out.evt.h5
```

where `files.txt` should be a textfile containing the full paths to the event HDF5 files to be merged, formatted with one file per line, `out.evt.h5` is the name of the merged output file, and `<N>` is the number of MPI ranks requested. For an HPC CPU node, `<N>` should be set to the number of CPU cores requested; in general, running on a personal computer is not recommended, and `<N>` should typically be <= 4.
