# -----------------------------------------------------------------
# PIHM
# Version: V 2.3
# Date: September 24, 2014 
# -----------------------------------------------------------------
# Programmer: Yuning Shi (yshi@psu.edu)
# -----------------------------------------------------------------
# Makefile for PIHM 
# -----------------------------------------------------------------

CC = gcc
CFLAGS = -g -O0

SRCDIR = ./src

SUNDIALS_PATH = /gpfs/home/yzs123/work/lib/sundials-2.2.0

SRCS_ =  pihm.c f.c read_alloc.c initialize.c update.c print.c is_sm_et.c \
	    f_function.c
HEADERS_ = pihm.h

EXECUTABLE = pihm
MSG = "...  Compiling PIHM  ..."

# Available modules are:
# FLUX-PIHM
# BGC	(Note: When BGC is declared, Flux-PIHM module will also be turned on
# 	because soil temperature needs to be simulated)

MODULE = BGC

ifeq ($(MODULE), FLUX-PIHM)
  SFLAGS = -D_FLUX_PIHM_ 
  MODULE_SRCS_ = noah/coupling.c noah/module_sf_noahlsm.c spa/spa.c \
		    noah/lsm_func.c
  MODULE_HEADERS_ =  noah/noah.h noah/flux_pihm.h spa/spa.h
  EXECUTABLE = flux-pihm
  MSG = "... Compiling FLUX-PIHM ..."
endif

ifeq ($(MODULE), BGC)
  SFLAGS = -D_BGC_ -D_FLUX_PIHM_ 
  MODULE_SRCS_ = time_func.c noah/coupling.c noah/module_sf_noahlsm.c \
		    spa/spa.c noah/lsm_func.c bgc/BGC_func.c \
		    bgc/presim_state_init.c bgc/make_zero_flux_struct.c \
		    bgc/restart_io.c bgc/firstday.c bgc/zero_srcsnk.c \
		    bgc/daily_bgc.c bgc/get_co2.c bgc/get_ndep.c \
		    bgc/precision_control.c bgc/daymet.c bgc/radtrans.c \
		    bgc/maint_resp.c bgc/phenology.c bgc/soilpsi.c \
		    bgc/daily_allocation.c bgc/canopy_et.c \
		    bgc/photosynthesis.c bgc/decomp.c bgc/annual_rates.c \
		    bgc/growth_resp.c bgc/state_update.c
  MODULE_HEADERS_ =  noah/noah.h noah/flux_pihm.h spa/spa.h bgc/bgc.h 
  EXECUTABLE = pihm-bgc
  MSG = "... Compiling PIHM-BGC ..."
endif

LIBS = -lm
INCLUDES = -I${SUNDIALS_PATH}/include -I${SUNDIALS_PATH}/include/cvode -I${SUNDIALS_PATH}/include/sundials
LFLAGS = -L${SUNDIALS_PATH}/lib -lsundials_cvode -lsundials_nvecserial

SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_))
MODULE_SRCS = $(patsubst %,$(SRCDIR)/%,$(MODULE_SRCS_))
HEADERS = $(patsubst %,$(SRCDIR)/%,$(HEADERS_))
MODULE_HEADERS = $(patsubst %,$(SRCDIR)/%,$(MODULE_HEADERS_))

#OBJS_ = $(SRCS_:.c=.o)
#OBJS = $(patsubst %,$(OBJDIR)/%,$(OBJS_))
#MODULE_OBJS_ = $(MODULE_SRCS_:.c=.o)
#MODULE_OBJS = $(patsubst %,$(OBJDIR)/%,$(MODULE_OBJS_))
OBJS = $(SRCS:.c=.o)
MODULE_OBJS = $(MODULE_SRCS:.c=.o)
  

all:	$(EXECUTABLE)
	@echo
	@echo $(MSG)
	@echo
#	@rm $(OBJDIR)/*.o
$(EXECUTABLE): $(OBJS) $(MODULE_OBJS)
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

#$(OBJDIR)/%.o: $(SRCDIR)/%.c
%.o: %.c
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -c $<  -o $@

tool:
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o convert src/tool/convert.c
#$(OBJS): | $(OBJDIR)
     
#$(OBJDIR):
#	mkdir $(OBJDIR)
.PHONY: clean

clean:
	$(RM) $(SRCDIR)/*.o $(SRCDIR)/noah/*.o $(SRCDIR)/spa/*.o *~ pihm flux-pihm
#	$(RM) $(OBJDIR)/*.o *~ pihm flux-pihm
