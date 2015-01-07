#-----------------------------------------------------------------
# PIHM Makefile
# Version: 2.4
# Date: Jan 10, 2014 
# -----------------------------------------------------------------
# Programmer: Yuning Shi (yshi@psu.edu)
# -----------------------------------------------------------------
# Makefile for PIHM 
# -----------------------------------------------------------------

CC = gcc
CFLAGS = -g -O0

SUNDIALS_PATH = /gpfs/home/yzs123/work/lib/sundials-2.2.0

SRCDIR = ./src
LIBS = -lm
INCLUDES = -I${SUNDIALS_PATH}/include -I${SUNDIALS_PATH}/include/cvode \
	    -I${SUNDIALS_PATH}/include/sundials
LFLAGS = -L${SUNDIALS_PATH}/lib -lsundials_cvode -lsundials_nvecserial

SRCS_ =  pihm.c f.c read_alloc.c initialize.c update.c print.c is_sm_et.c \
	    f_function.c forcing.c
HEADERS_ = pihm.h
EXECUTABLE = pihm
MSG = "...  Compiling PIHM  ..."

FLUX_SFLAGS = -D_FLUX_PIHM_ 
FLUX_SRCS_ = noah/coupling.c noah/module_sf_noahlsm.c spa/spa.c \
		noah/lsm_func.c
FLUX_HEADERS_ =  noah/noah.h spa/spa.h
FLUX_EXECUTABLE = flux-pihm
FLUX_MSG = "... Compiling FLUX-PIHM ..."

BGC_SFLAGS = -D_BGC_ -D_FLUX_PIHM_ 
BGC_SRCS_ = time_func.c noah/coupling.c noah/module_sf_noahlsm.c spa/spa.c \
		    noah/lsm_func.c bgc/BGC_func.c bgc/presim_state_init.c \
		    bgc/make_zero_flux_struct.c bgc/restart_io.c \
		    bgc/firstday.c bgc/zero_srcsnk.c bgc/daily_bgc.c \
		    bgc/get_co2.c bgc/get_ndep.c bgc/precision_control.c \
		    bgc/daymet.c bgc/radtrans.c bgc/maint_resp.c \
		    bgc/phenology.c bgc/soilpsi.c bgc/daily_allocation.c \
		    bgc/canopy_et.c bgc/photosynthesis.c bgc/decomp.c \
		    bgc/annual_rates.c bgc/growth_resp.c bgc/state_update.c \
		    bgc/mortality.c bgc/check_balance.c bgc/summary.c \
		    bgc/metarr_init.c bgc/bgc_spinup.c
BGC_HEADERS_ =  noah/noah.h noah/flux_pihm.h spa/spa.h bgc/bgc.h 
BGC_EXECUTABLE = pihm-bgc
BGC_MSG = "... Compiling PIHM-BGC ..."

SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_))
HEADERS = $(patsubst %,$(SRCDIR)/%,$(HEADERS_))
OBJS = $(SRCS:.c=.o)

FLUX_SRCS = $(patsubst %,$(SRCDIR)/%,$(FLUX_SRCS_))
FLUX_HEADERS = $(patsubst %,$(SRCDIR)/%,$(FLUX_HEADERS_))
FLUX_OBJS = $(FLUX_SRCS:.c=.o)

BGC_SRCS = $(patsubst %,$(SRCDIR)/%,$(BGC_SRCS_))
BGC_HEADERS = $(patsubst %,$(SRCDIR)/%,$(BGC_HEADERS_))
BGC_OBJS = $(BGC_SRCS:.c=.o)

pihm: $(OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(LFLAGS) $(LIBS)

flux-pihm: $(OBJS) $(FLUX_OBJS)
	@echo
	@echo $(FLUX_MSG)
	@echo
	@$(CC) $(CFLAGS) $(FLUX_SFLAGS) $(INCLUDES) -o $(FLUX_EXECUTABLE) $(OBJS) $(FLUX_OBJS) $(LFLAGS) $(LIBS)
tool:
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o convert src/tool/convert.c

$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(FLUX_OBJS): %.o: %.c
	$(CC) $(CFLAGS) $(FLUX_SFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	@echo
	@echo "... Cleaning ..."
	@echo
	@$(RM) $(SRCDIR)/*.o $(SRCDIR)/*/*.o *~ $(EXECUTABLE) $(FLUX_EXECUTABLE) $(BGC_EXECUTABLE)
#	$(RM) $(OBJDIR)/*.o *~ pihm flux-pihm

