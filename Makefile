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

SRCS_ =  pihm.c f.c read_alloc.c initialize.c update.c print.c is_sm_et.c f_function.c
HEADERS_ = pihm.h

EXECUTABLE = pihm
MSG = "...  Compiling PIHM  ..."

# Available modules are:
# FLUX-PIHM
# BBGC

MODULE = BBGC

ifeq ($(MODULE), FLUX-PIHM)
  SFLAGS = -D_FLUX_PIHM_ 
  MODULE_SRCS_ = noah/coupling.c noah/module_sf_noahlsm.c spa/spa.c noah/lsm_func.c
  MODULE_HEADERS_ =  noah/noah.h noah/flux_pihm.h spa/spa.h
  EXECUTABLE = flux-pihm
  MSG = "... Compiling FLUX-PIHM ..."
endif

ifeq ($(MODULE), BBGC)
  SFLAGS = -D_BBGC_ 
  MODULE_SRCS_ = bgc/coupling.c 
  MODULE_HEADERS_ =  bgc/pihm_bgc.h 
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

#$(OBJS): | $(OBJDIR)
     
#$(OBJDIR):
#	mkdir $(OBJDIR)
.PHONY: clean

clean:
	$(RM) $(SRCDIR)/*.o $(SRCDIR)/noah/*.o $(SRCDIR)/spa/*.o *~ pihm flux-pihm
#	$(RM) $(OBJDIR)/*.o *~ pihm flux-pihm
