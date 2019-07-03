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

#define VERSION    "0.13.0-alpha"

/*
 * SUNDIAL Header Files
 */
/* CVODE header file */
#include "cvode.h"

/* CVSPGMR linear header file */
#include "cvode_spgmr.h"

/* Definition of type N_Vector */
#if defined(_CVODE_OMP)
# include "nvector_openmp.h"
#else
# include "nvector_serial.h"
#endif

/* UnitRoundoff, RSqrt, SQR functions */
#include "sundials_math.h"

/* CVDENSE header file */
#include "cvode_dense.h"

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
