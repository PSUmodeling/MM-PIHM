#-----------------------------------------------------------------
# MM-PIHM Makefile
# -----------------------------------------------------------------

CC = gcc
CFLAGS = -g -O2

ifeq ($(WARNING), on)
  CFLAGS += -Wall -Wextra
endif

ifeq ($(DEBUG), on)
  CFLAGS += -O0
endif

ifneq ($(OMP), off)
  CFLAGS += -fopenmp
endif

CMAKETEST=$(shell cmake --version 2> /dev/null)

ifeq ($(CMAKETEST),)
  CMAKE_EXIST = 0
  OS := $(shell uname)
ifeq ($(OS),Darwin)
    CMAKE_VERS = cmake-3.7.2-Darwin-x86_64
    CMAKE = $(PWD)/$(CMAKE_VERS)/CMake.app/Contents/bin/cmake
else
    CMAKE_VERS = cmake-3.7.2-Linux-x86_64
    CMAKE = $(PWD)/$(CMAKE_VERS)/bin/cmake
endif
else
  CMAKE_EXIST = 1
  CMAKE=cmake
endif

CVODE_PATH = ./cvode/instdir

SRCDIR = ./src
LIBS = -lm -Wl,-rpath,$(CVODE_PATH)/lib
INCLUDES = \
	-I$(SRCDIR)/include\
	-I$(CVODE_PATH)/include\
	-I$(CVODE_PATH)/include/cvode\
	-I$(CVODE_PATH)/include/sundials\
	-I$(CVODE_PATH)/include/nvector

LFLAGS = -lsundials_cvode -L$(CVODE_PATH)/lib
ifeq ($(CVODE_OMP), on)
  LFLAGS += -lsundials_nvecopenmp
else
  LFLAGS += -lsundials_nvecserial
endif

SFLAGS = -D_PIHM_

ifeq ($(CVODE_OMP), on)
  SFLAGS += -D_CVODE_OMP
endif

ifeq ($(DEBUG), on)
  SFLAGS += -D_DEBUG_
endif

SRCS_ = main.c\
	custom_io.c\
	forcing.c\
	free_mem.c\
	hydrol.c\
	init_forc.c\
	init_lc.c\
	init_mesh.c\
	init_river.c\
	init_soil.c\
	init_topo.c\
	initialize.c\
	is_sm_et.c\
	lat_flow.c\
	map_output.c\
	ode.c\
	optparse.c\
	pihm.c\
	print.c\
	read_alloc.c\
	read_att.c\
	read_bc.c\
	read_calib.c\
	read_forc.c\
	read_func.c\
	read_ic.c\
	read_lai.c\
	read_lc.c\
	read_mesh.c\
	read_para.c\
	read_river.c\
	read_soil.c\
	read_tecplot.c\
	river_flow.c\
	soil.c\
	spinup.c\
	time_func.c\
	update.c\
	util_func.c\
	vert_flow.c

HEADERS_ = \
	include/elem_struct.h\
	include/pihm_const.h\
	include/pihm_func.h\
	include/pihm_input_struct.h\
	include/pihm_struct.h\
	include/pihm.h\
	include/river_struct.h

MODULE_HEADERS_ =
EXECUTABLE = pihm
MSG = "...  Compiling PIHM  ..."

#-------------------
# PIHM-FBR
#-------------------
ifeq ($(MAKECMDGOALS),pihm-fbr)
  SFLAGS += -D_FBR_
  MODULE_SRCS_ = \
	fbr/init_geol.c\
	fbr/read_bedrock.c\
  	fbr/read_geol.c
  MODULE_HEADERS_ =
  EXECUTABLE = pihm-fbr
  MSG = "... Compiling PIHM-FBR ..."
endif

#-------------------
# Flux-PIHM
#-------------------
ifeq ($(MAKECMDGOALS),flux-pihm)
  SFLAGS += -D_NOAH_
  MODULE_SRCS_ = \
	noah/lsm_init.c\
	noah/lsm_func.c\
	noah/lsm_read.c\
	noah/noah.c\
	noah/noah_glacial_only.c\
	noah/topo_radn.c\
	spa/spa.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm
  MSG = "... Compiling Flux-PIHM ..."
endif

#-------------------
# Flux-PIHM-FBR
#-------------------
ifeq ($(MAKECMDGOALS),flux-pihm-fbr)
  SFLAGS += -D_NOAH_ -D_FBR_
  MODULE_SRCS_ = \
  	fbr/init_geol.c\
	fbr/read_bedrock.c\
	fbr/read_geol.c\
	noah/lsm_init.c\
	noah/noah.c\
	noah/noah_glacial_only.c\
	noah/lsm_func.c\
	noah/lsm_read.c\
	noah/topo_radn.c\
	spa/spa.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm-fbr
  MSG = "... Compiling Flux-PIHM-FBR ..."
endif

#-------------------
# RT-Flux-PIHM
#-------------------
ifeq ($(MAKECMDGOALS),rt-flux-pihm)
  SFLAGS += -D_RT_ -D_NOAH_
  MODULE_SRCS_=\
	noah/lsm_init.c\
	noah/lsm_func.c\
	noah/lsm_read.c\
	noah/noah.c\
	noah/noah_glacial_only.c\
	noah/topo_radn.c\
	spa/spa.c\
	rt/read_chem.c\
	rt/read_prep.c\
	rt/rt_util.c
	#rt/flux_trans.c\
	#rt/lookup.c\
	#rt/os3d.c\
	#rt/react.c\
	#rt/read_cini.c\
	#rt/rt.c\
	#rt/speciation.c
  MODULE_HEADERS_ =\
	include/spa.h
  EXECUTABLE = rt-flux-pihm
  MSG = "... Compiling RT-Flux-PIHM ..."
endif

#-------------------
# Flux-PIHM-BGC
#-------------------
ifeq ($(MAKECMDGOALS),flux-pihm-bgc)
  SFLAGS += -D_NOAH_ -D_BGC_ -D_DAILY_
ifeq ($(LUMPED), on)
  SFLAGS += -D_LUMPED_
endif
ifeq ($(LEACHING), on)
  SFLAGS += -D_LEACHING_
endif
  MODULE_SRCS_= \
	bgc/bgc_init.c\
	bgc/bgc_read.c\
	bgc/bgc.c\
	bgc/canopy_cond.c\
	bgc/check_balance.c\
	bgc/daily_allocation.c\
	bgc/decomp.c\
	bgc/firstday.c\
	bgc/get_co2.c\
	bgc/get_ndep.c\
	bgc/growth_resp.c\
	bgc/maint_resp.c\
	bgc/make_zero_flux_struct.c\
	bgc/mortality.c\
	bgc/ntransport.c\
	bgc/phenology.c\
	bgc/photosynthesis.c\
	bgc/precision_control.c\
	bgc/presim_state_init.c\
	bgc/radtrans.c\
	bgc/restart_io.c\
	bgc/soilpsi.c\
	bgc/state_update.c\
	bgc/summary.c\
	bgc/zero_srcsnk.c\
	noah/daily.c\
	noah/lsm_func.c\
	noah/lsm_init.c\
	noah/lsm_read.c\
	noah/noah.c\
	noah/noah_glacial_only.c\
	noah/topo_radn.c\
	spa/spa.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm-bgc
  MSG = "... Compiling Flux-PIHM-BGC ..."
endif

#-------------------
# Flux-PIHM-Cycles
#-------------------
CYCLES_PATH = ../Cycles_esm/src
RQD_CYCLES_VERS = R0.8.0-alpha
ifeq ($(MAKECMDGOALS),flux-pihm-cycles)
  SFLAGS += -D_NOAH_ -D_CYCLES_ -D_DAILY_
  MODULE_SRCS_= \
	cycles/cycles.c\
  	cycles/cycles_read.c\
	cycles/cycles_init.c\
	cycles/ntransport.c\
	cycles/update_prof.c\
	noah/daily.c\
	noah/lsm_func.c\
	noah/lsm_init.c\
	noah/lsm_read.c\
	noah/noah.c\
	noah/noah_glacial_only.c\
	noah/topo_radn.c\
	spa/spa.c
  CYCLES_SRCS_ =\
	crop.c\
	crop_harvest.c\
	crop_process.c\
	crop_thermal_time.c\
	crop_transpiration.c\
	daily_operation.c\
	fertilization.c\
	field_operation.c\
	growing_crop.c\
	irrigation.c\
	make_zero_flux_struct.c\
	read_crop.c\
	read_operation.c\
	residue.c\
	restart.c\
	soil_carbon.c\
	soil_nitrogen.c\
	soil_solute.c\
	tillage.c\
	time_func.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm-cycles
  MSG = "... Compiling Flux-PIHM-Cycles ..."
endif

SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_))
HEADERS = $(patsubst %,$(SRCDIR)/%,$(HEADERS_))
OBJS = $(SRCS:.c=.o)

MODULE_SRCS = $(patsubst %,$(SRCDIR)/%,$(MODULE_SRCS_))
MODULE_HEADERS = $(patsubst %,$(SRCDIR)/%,$(MODULE_HEADERS_))
MODULE_OBJS = $(MODULE_SRCS:.c=.o)

CYCLES_SRCS = $(patsubst %,$(CYCLES_PATH)/%,$(CYCLES_SRCS_))
CYCLES_OBJS = $(CYCLES_SRCS:.c=.o)

.PHONY: all clean help cvode cmake

help:			## Show this help
	@echo
	@echo "Makefile for MM-PIHM"
	@echo
	@echo "USAGE:"
	@echo
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'
	@echo
	@echo "NOTE: Please always \"make clean\" when switching from one module to another!"
	@echo

all:			## Install cvode and compile PIHM
all:	cvode pihm

cmake:
ifneq ($(CMAKE_EXIST),1)
	@echo "Download CMake from cmake.org"
	@curl https://cmake.org/files/v3.7/$(CMAKE_VERS).tar.gz -o $(CMAKE_VERS).tar.gz &> /dev/null
	@echo
	@echo "Extract $(CMAKE_VERS).tar.gz"
	@tar xzf $(CMAKE_VERS).tar.gz
endif

cvode:			## Install cvode library
cvode:	cmake
	@echo "Install CVODE library"
	@cd cvode && mkdir -p instdir && mkdir -p builddir
	@cd $(CVODE_PATH) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=../instdir -DEXAMPLES_ENABLE=OFF -DEXAMPLES_INSTALL=OFF ../
	@cd $(CVODE_PATH) && make && make install
	@echo "CVODE library installed."
ifneq ($(CMAKE_EXIST),1)
	@echo "Remove CMake files"
	@$(RM) -r $(CMAKE_VERS).tar.gz $(CMAKE_VERS)
endif

pihm:			## Compile PIHM
pihm:	$(OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(LFLAGS) $(LIBS)

pihm-fbr:		## Compile PIHM-FBR (PIHM with fractured bedrock module)
pihm-fbr: $(OBJS) $(MODULE_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

flux-pihm:		## Compile Flux-PIHM (PIHM with land surface module, adapted from Noah LSM)
flux-pihm: $(OBJS) $(MODULE_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

flux-pihm-fbr:		## Compile Flux-PIHM-FBR (PIHM with land surface and fractured bedrock modules)
flux-pihm-fbr: $(OBJS) $(MODULE_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

rt-flux-pihm:		## Compile RT-Flux-PIHM (PIHM with land surface and reactive transport modules)
rt-flux-pihm: $(OBJS) $(MODULE_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

flux-pihm-bgc:		## Compile Flux-PIHM-BGC (Flux-PIHM with Biogeochemical module, adapted from Biome-BGC)
flux-pihm-bgc: $(OBJS) $(MODULE_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

flux-pihm-cycles:	## Compile PIHM-Cycles (Flux-PIHM with crop module, adapted from Cycles)
flux-pihm-cycles: check_cycles_vers $(OBJS) $(MODULE_OBJS) $(CYCLES_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(CYCLES_OBJS) $(LFLAGS) $(LIBS)

check_cycles_vers:
	@util/check_cycles_vers.sh $(CYCLES_PATH) $(RQD_CYCLES_VERS)

%.o: %.c $(HEADERS) $(MODULE_HEADERS)
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -c $<  -o $@


clean:			## Clean executables and objects
	@echo
	@echo "... Cleaning ..."
	@echo
	@$(RM) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(CYCLES_PATH)/*.o *~ pihm pihm-fbr flux-pihm flux-pihm-fbr rt-flux-pihm flux-pihm-bgc flux-pihm-cycles rt-flux-pihm
