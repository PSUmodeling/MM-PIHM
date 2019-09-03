#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <stdarg.h>
#include <sys/stat.h>
#if defined(_WIN32) || defined(_WIN64)
# include <windows.h>
# include <direct.h>
# include <io.h>
#else
# include <unistd.h>
#endif
#if defined(_OPENMP)
# include <omp.h>
#endif

#define VERSION    "0.14.0-alpha"

/*
 * SUNDIAL Header Files
 */
/* Prototypes for CVODE fcts., consts. */
#include "cvode/cvode.h"

/* Access to SPGMR SUNLinearSolver */
#include "sunlinsol/sunlinsol_spgmr.h"

/* Access to N_Vector */
#if defined(_CVODE_OMP)
# include "nvector/nvector_openmp.h"
#else
# include "nvector/nvector_serial.h"
#endif

/* Definition of macros SUNSQR and EXP */
#include "sundials/sundials_math.h"

/* Prototypes for small dense fcts. */
#include "sundials/sundials_dense.h"

#if defined(_NOAH_)
# include "spa.h"
#endif

#include "custom_io.h"

#include "pihm_const.h"
#include "pihm_input_struct.h"
#include "elem_struct.h"
#include "river_struct.h"
#include "pihm_struct.h"
#include "pihm_func.h"

#endif
