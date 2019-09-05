MM-PIHM [![Build Status](https://travis-ci.org/PSUmodeling/MM-PIHM.svg?branch=master)](https://travis-ci.org/PSUmodeling/MM-PIHM)
=======

The Multi-Modular Penn State Integrated Hydrologic Model (MM-PIHM) is a physically based watershed model with multiple optional modules.
MM-PIHM is the **sweetest** PIHM, ever!

The current release contains the source code for PIHM, Flux-PIHM, Flux-PIHM-BGC, and RT-Flux-PIHM.

PIHM is a spatially-distributed, physically based hydrologic model.
Flux-PIHM adds a land surface model (adapted from the Noah land surface model) to PIHM for
the simulation of land surface processes.
Flux-PIHM-BGC couples Flux-PIHM with a terrestrial ecosystem model (adapted from Biome-BGC) that enables the simulation of carbon and nitrogen cycles.
RT-Flux-PIHM couples Flux-PIHM with a multicomponent subsurface reactive transport module.
A deep groundwater module can be turned on for PIHM, Flux-PIHM, and RT-Flux-PIHM.

MM-PIHM is open source software licensed under the MIT License.
All bug reports and feature requests should be submitted using the [Issues](https://github.com/PSUmodeling/MM-PIHM/issues) page.

## Usage

The following guide applies to UNIX (include Mac OS) systems.
For instructions on how to install MM-PIHM on Windows, please refer to this [guide](https://gist.github.com/shiyuning/867d5af0a3a6345b50ec1b193a71e4be).

### Installing CVODE

MM-PIHM uses the SUNDIALS CVODE v2.9.0 implicit solvers.
The CVODE Version 2.9.0 source code is provided with the MM-PIHM package for users' convenience.
SUNDIALS (:copyright:2012--2016) is copyrighted software produced at the Lawrence Livermore National Laboratory.
A SUNDIALS copyright note can be found in the `cvode` directory.

If you already have CVODE v2.9.0 installed, you can edit the Makefile and point `CVODE_PATH` to your CVODE directory.
Otherwise, you need to install CVODE before compiling MM-PIHM, by doing

```shell
$ make cvode
```

in your MM-PIHM directory.

Currently CMake (version 2.8.1 or higher) is the only supported method of CVODE installation.
If CMake is not available on your system, the CMake Version 3.7.2 binary for Linux (or Mac OS, depending on your OS) will be downloaded from [http://www.cmake.org](http://www.cmake.org) automatically when you choose to `make cvode`.

### Installing MM-PIHM

Once CVODE is installed, you can compile MM-PIHM models from the MM-PIHM directory by doing

```shell
$ make [model]
```

The `[model]` should be replaced by the name of model that you want to compile, which could be `pihm`, `flux-pihm`, `flux-pihm-bgc`, or `rt-flux-pihm`.

The command

```shell
$ make clean
```

will clean the executables and object files.

Note: If you want to switch from one MM-PIHM model to another one, you must `make clean` first.

A help message will appear if you run `make`.

#### Installation options

The deep groundwater module (DGW) can be turned on during compilation.
To turn on the DGW module, compile MM-PIHM models using

```shell
$ make DGW=on [model]
```

By default, MM-PIHM is paralleled using OpenMP, which significantly improves the computational efficiency of MM-PIHM models.
CVODE, however, is not implemented using OpenMP by default.
According to CVODE document, CVODE state variables (i.e., unknowns) "should be of length at least 100, 000 before the overhead associated with creating and using the threads is made up by the parallelism in the vector calculations".
In other words, you should using OpenMP for CVODE if your model domain has about 30, 000 or more model grids.
If you do want to test using OpenMP for CVODE, you can compile MM-PIHM models using

```shell
$ make CVODE_OMP=on [model]
```

Note that in order to use OpenMP for CVODE, you also need to turn on the OPENMP_ENABLE option when using CMake to install CVODE.

You can also turn off OpenMP for MM-PIHM (**NOT RECOMMENDED**):

```shell
$ make OMP=off [model]
```

By default, PIHM compilation is optimized using the `-O0` gcc option.
If you wish to further accelerate the model (**NOT RECOMMENDED**), you may want to use

```shell
$ make DEBUG=off [model]
```

which will compile using `-O2` gcc option.

### Running MM-PIHM

#### Setting up OpenMP environment

To optimize PIHM efficiency, you need to set the number of threads in OpenMP.
For example, in command line

```shell
$ export OMP_NUM_THREADS=12
```

The command above will enable MM-PIHM model simulations using twelve (12) OpenMP threads.

If you use a PBS script, you must require the right number of ppn (processor cores per node) before setting the number of threads.
The ppn should be the same as the number of threads you want to use.
For example, your PBS script may look like

```shell
#PBS -l nodes=1:ppn=12
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -l pmem=1gb

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=12
./pihm example
```

#### Running MM-PIHM models

Now you can run MM-PIHM models using:

```shell
$ ./[model] [-b] [-c] [-d] [-f] [-s] [-t] [-V] [-v] [-o dir_name] [project]
```

where `[model]` is the installed executable, `[project]` is the name of the project, and `[-bcdfostVv]` are optional parameters.

The optional `-b` parameter will turn on the brief mode with minimum screen output.

The optional `-c` parameter will turn on the elevation correction mode.
Surface elevation of all model grids will be checked, and changed if needed before simulation, to avoid surface sinks.

The optional `-d` parameter will turn on the debug mode.
In debug mode, helpful information is displayed on screen and a CVODE log file will be produced.

The optional `-f` parameter will turn on fixed length spin-up mode, in which spin-up simulations will only stop at the specified maximum spin-up years, but not at equilibrium.

The optional `-s` parameter will turn on silence mode without screen output during simulations.

The optional `-t` parameter will turn on Tecplot output.

The `-V` parameter will display model version.
Note that model will quit after displaying the version information.
No simulation will be performed when using the `-V` parameter.

The optional `-v` parameter will turn on the verbose mode.

The optional `-o` parameter will specify the name of directory to store model output.
All model output variables will be stored in the `output/dir_name` directory when `-o` option is used.
If `-o` parameter is not used, model output will be stored in a directory named after the project and the system time when the simulation is executed.

Example input files are provided with each release.
For a description of input files, please refer to the *User's Guide*
 that can be downloaded from the [release page](https://github.com/PSUmodeling/MM-PIHM/releases).

### Penn State Users

The Penn State ICS ACI system support both batch job submissions and interactive jobs.
The clusters for batch job submissions, ACI-B, usually support twelve processors per node for OMP jobs.
The cluster for interactive jobs, ACI-I, supports twenty processors per node for OMP jobs, but is limited to short jobs only.
