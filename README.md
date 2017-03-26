MM-PIHM [![Build Status](https://travis-ci.org/PSUmodeling/MM-PIHM.svg?branch=master)](https://travis-ci.org/PSUmodeling/MM-PIHM)
=======

The Multi-Modular Penn State Integrated Hydrologic Model is a physics based hydrologic model with multiple optional modules, including land surface, reactive transport, and biogeochemistry modules.
MM-PIHM is the **sweetest** PIHM, ever!

The first release includes PIHM, Flux-PIHM, Flux-PIHM-BGC, and Flux-PIHM-EnKF.
Future release will include RT-Flux-PIHM.

## Usage

### Install Cvode

MM-PIHM uses the SUNDIALS  CVODE v2.9.0 implicit solvers.
If you already have CVODE v2.9.0 installed, you can edit the Makefile and point `CVODE_PATH` to your CVODE directory.
Otherwise, you need to install CVODE before compiling MM-PIHM.
Currently CMake (version 2.8.1 or higher) is the only supported method of cvode installation.
You can download the latest version of CMake from [http://www.cmake.org](http://www.cmake.org), and install.
Note: you may need to add Cmake directory to your `PATH` variable, by adding

```shell
PATH=$PATH:YOUR_CMAKE_DIRECTORY
```

to your `.bash_profile` file in your home directory, and then `source ~/.bash_profile`.

Once you have Cmake installed, you can install CVODE using

```shell
$ make cvode
```

in your MM-PIHM directory.
