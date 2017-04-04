#-----------------------------------------------------------------
# MM-PIHM Makefile
# -----------------------------------------------------------------

CC = gcc
CFLAGS = -g -O2

ifeq ($(DEBUG), on)
CFLAGS += -O0 -Wall -Wextra
endif

ifneq ($(OMP), off)
CFLAGS += -fopenmp
endif

CVODE_PATH = ./cvode/instdir

SRCDIR = ./src
LIBS = -lm -Wl,-rpath,${CVODE_PATH}/lib
INCLUDES = \
	-I${SRCDIR}/include\
	-I${CVODE_PATH}/include\
	-I${CVODE_PATH}/include/cvode\
	-I${CVODE_PATH}/include/sundials\
	-I${CVODE_PATH}/include/nvector


LFLAGS = -lsundials_cvode -L${CVODE_PATH}/lib
ifeq ($(CVODE_OMP), on)
LFLAGS += -lsundials_nvecopenmp
else
LFLAGS += -lsundials_nvecserial
endif

SFLAGS = -D_PIHM_
ifeq ($(CVODE_OMP), on)
SFLAGS += -D_CVODE_OMP
endif

SRCS_ = main.c\
	forcing.c\
	hydrol.c\
	initialize.c\
	is_sm_et.c\
	lat_flow.c\
	misc_func.c\
	pihm.c\
	print.c\
	read_alloc.c\
	read_func.c\
	river_flow.c\
	soil.c\
	update.c\
	vert_flow.c

HEADERS_ = \
	include/elem_struct.h\
	include/pihm.h\
	include/pihm_const.h\
	include/pihm_func.h\
	include/pihm_input_struct.h\
	include/pihm_struct.h\
	include/river_struct.h

MODULE_HEADERS_ =
EXECUTABLE = pihm
MSG = "...  Compiling PIHM  ..."

#-------------------
# Flux-PIHM
#-------------------
ifeq ($(MAKECMDGOALS),flux-pihm)
  SFLAGS = -D_PIHM_ -D_NOAH_ 
  MODULE_SRCS_ = \
  	noah/lsm_func.c\
	noah/lsm_init.c\
  	noah/lsm_read.c\
	noah/noah.c\
	spa/spa.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm
  MSG = "... Compiling Flux-PIHM ..."
endif

#-------------------
# RT-Flux-PIHM
#-------------------
#ifeq ($(MAKECMDGOALS),rt-flux-pihm)
#  SFLAGS = -D_PIHM_ -D_RT_ -D_NOAH_
#  MODULE_SRCS_=\
#  	noah/coupling.c\
#	noah/module_sf_noahlsm.c\
#	spa/spa.c\
#	noah/lsm_func.c\
#	rt/rt.c\
#	rt/react.c\
#	rt/os3d.c
#  MODULE_HEADERS_ =\
#	spa/spa.h\
#	rt/rt.h
#  EXECUTABLE = rt-flux-pihm
#  MSG = "... Compiling RT-Flux-PIHM ..."
#endif

#-------------------
# Flux-PIHM-BGC
#-------------------
ifeq ($(MAKECMDGOALS),flux-pihm-bgc)
  SFLAGS = -D_PIHM_ -D_NOAH_ -D_BGC_ -D_DAILY_
  MODULE_SRCS_= \
	bgc/annual_rates.c\
	bgc/bgc_init.c\
	bgc/bgc_read.c\
	bgc/bgc_spinup.c\
	bgc/canopy_cond.c\
	bgc/check_balance.c\
	bgc/daily_allocation.c\
	bgc/daily_bgc.c\
	bgc/daymet.c\
	bgc/decomp.c\
	bgc/firstday.c\
	bgc/get_co2.c\
	bgc/get_ndep.c\
	bgc/growth_resp.c\
	bgc/maint_resp.c\
	bgc/make_zero_flux_struct.c\
	bgc/metarr_init.c\
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
	spa/spa.c
  MODULE_HEADERS_ = include/spa.h
  EXECUTABLE = flux-pihm-bgc
  MSG = "... Compiling Flux-PIHM-BGC ..."
endif

#-------------------
# Flux-PIHM-Cycles
#-------------------
CYCLES_PATH = ../Cycles/src
ifeq ($(MAKECMDGOALS),flux-pihm-cycles)
  SFLAGS = -D_PIHM_ -D_NOAH_ -D_CYCLES_ -D_DAILY_
  MODULE_SRCS_= \
	cycles/cycles_func.c\
	cycles/cycles_init.c\
  	cycles/cycles_read.c\
	noah/daily.c\
	noah/lsm_func.c\
	noah/lsm_init.c\
	noah/lsm_read.c\
  	noah/noah.c\
	spa/spa.c
  CYCLES_SRCS_ = \
	Crop.c\
	CropHarvest.c\
	CropProcess.c\
	CropThermalTime.c\
	CropTranspiration.c\
	DailyOperation.c\
	Fertilization.c\
	FieldOperation.c\
	Irrigation.c\
	Residue.c\
	Soil.c\
	SoilCarbon.c\
	SoilEvaporation.c\
	SoilNitrogen.c\
	SoilSolute.c\
	Tillage.c
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

.PHONY: all clean help cvode

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

cvode:			## Install cvode library
cvode:
	@cd cvode && mkdir -p instdir && mkdir -p builddir
	@cd ${CVODE_PATH} && cmake -DCMAKE_INSTALL_PREFIX=../instdir -DEXAMPLES_ENABLE=OFF -DEXAMPLES_INSTALL=OFF ../
	@cd ${CVODE_PATH} && make && make install
	@echo "SUNDIALS library installed."

pihm:			## Compile PIHM
pihm:	$(OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(LFLAGS) $(LIBS)

flux-pihm:		## Complile Flux-PIHM (PIHM with land surface module, adapted from Noah LSM)
flux-pihm: $(OBJS) $(MODULE_OBJS)
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
flux-pihm-cycles: $(OBJS) $(MODULE_OBJS) $(CYCLES_OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(CYCLES_OBJS) $(LFLAGS) $(LIBS)

%.o: %.c $(HEADERS) $(MODULE_HEADERS)
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -c $<  -o $@


clean:			## Clean executables and objects
	@echo
	@echo "... Cleaning ..."
	@echo
	@$(RM) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(CYCLES_PATH)/*.o *~ pihm flux-pihm flux-pihm-bgc flux-pihm-cycles rt-flux-pihm
