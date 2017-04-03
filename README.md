MM-PIHM [![Build Status](https://travis-ci.org/PSUmodeling/MM-PIHM.svg?branch=master)](https://travis-ci.org/PSUmodeling/MM-PIHM)
=======

The Multi-Modular Penn State Integrated Hydrologic Model is a physics based hydrologic model with multiple optional modules, including land surface, reactive transport, and biogeochemistry modules.
MM-PIHM is the **sweetest** PIHM, ever!

The first release includes PIHM, Flux-PIHM, and Flux-PIHM-BGC.
Future release will include RT-Flux-PIHM.

## Usage

### Install CVODE

MM-PIHM uses the SUNDIALS CVODE v2.9.0 implicit solvers.
If you already have CVODE v2.9.0 installed, you can edit the Makefile and point `CVODE_PATH` to your CVODE directory.
Otherwise, you need to install CVODE before compiling MM-PIHM.
Currently CMake (version 2.8.1 or higher) is the only supported method of cvode installation.
You can download the latest version of CMake from [http://www.cmake.org](http://www.cmake.org), and install CMake.
Note: you may need to add CMake directory to your `PATH` variable, by adding

```shell
PATH=$PATH:YOUR_CMAKE_DIRECTORY
```

to your `.bash_profile` file in your home directory, and then `source ~/.bash_profile`.

Once you have Cmake installed, you can install CVODE using

```shell
$ make cvode
```

in your MM-PIHM directory.

### Install MM-PIHM

Once CVODE is installed, you can compile MM-PIHM models from the MM-PIHM directory by doing

```shell
$ make [model]
```

The `[model]` should be replaced by the name of model that you want to compile, which could be `pihm`, `flux-pihm`, or `flux-pihm-bgc`.
Note: If you want to switch from one MM-PIHM family model to another one, you always need to `make clean` before making another model.

#### Installation options

By default, PIHM is paralledled using OpenMP, which significantly improves the computational efficiency of MM-PIHM models, especially Flux-PIHM and Flux-PIHM-BGC.
CVODE, however, is not implemented using OpenMP by default.
According to CVODE document, CVODE state variables (unknowns) "should be of length at least 100, 000 before the overhead associated with creating and using the threads is made up by the parallelism in the vector calculations".
If you do want to test using OpenMP for CVODE, you can compile MM-PIHM models uing 

```shell
$ make CVODE_OMP=on [model]
```

You can also trun off OpenMP (NOT RECOMMENDED):

```shell
$ make OMP=off [model]
```

By default, PIHM compilation is optimized using `-O2` gcc option.
If you wish to debug using gdb, you may want to use

```shell
$ make DEBUG=on [model]
```

which will compile using `-O0` gcc option.

### Run MM-PIHM

#### Set up OpenMP environment

To optimize PIHM efficiency, you need to set the number of threads in OpenMP.
For example, in command line

```shell
$ export OMP_NUM_THREADS=8
```

The command above will enable MM-PIHM model simulaions using eight (8) OpenMP threads.

If you use a PBS script, you must require the right number of ppn (processor cores per node) before setting the number of threads.
The ppn should be the same as the number of threads you want to use.
For example, your PBS script may look like

```shell
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -l pmem=1gb

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
./pihm example
```

Now you can run MM-PIHM models:

```shell
$ ./<model> [-v] [-d] [-c] [-o dir_name] <project>
```

where `<model>` is the name of the MM-PIHM model, `<project>` is the name of the project, and [-vdco] are
optional parameters.

The optional `-v` parameter will turn on the verbose mode.
The optional `-d` parameter will turn on the debug mode.
The optional `-c` parameter will examine the surface elevation of all model elements and fix potential sinks.
The optional `-o` parameter will specify the name of directory to store model output.
All model output variables will be stored in the `output/dir_name` directory when `-o` option is used.
Otherwise, model output will be stored in a directory named after the project and the system time when the simulation is executed.
