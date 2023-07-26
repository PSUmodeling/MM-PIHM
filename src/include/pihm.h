#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <stdarg.h>
#include <float.h>
#include <sys/stat.h>
#if defined(_WIN32) || defined(_WIN64)
# include <windows.h>
# include <direct.h>
# include <io.h>
#else
# include <unistd.h>
#endif
#if defined(unix) || defined(__unix__) || defined(__unix)
# include <fenv.h>
#endif
#if defined(_OPENMP)
# include <omp.h>
#endif

#define VERSION             "1.0.0.post"

// SUNDIAL Header Files
#include "cvode/cvode.h"    // Prototypes for CVODE fcts., consts.
#include "sunlinsol/sunlinsol_spgmr.h"  // Access to SPGMR SUNLinearSolver
#if defined(_CVODE_OMP)
# include "nvector/nvector_openmp.h"    // Access to N_Vector
#else
# include "nvector/nvector_serial.h"
#endif
#include "sundials/sundials_math.h"     // Definition of macros SUNSQR and EXP
#include "sundials/sundials_dense.h"    // Prototypes for small dense fcts.

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
#include "pihm_errors.h"

#endif
