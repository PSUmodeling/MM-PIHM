#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>

#if !defined(_WIN32)
 #include <unistd.h>
#else
 #include <io.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* SUNDIAL Header Files */
#include "sundials_types.h"     /* realtype, integertype, booleantype
                                 * definition */
#include "cvode.h"              /* CVODE header file */
#include "cvode_spgmr.h"        /* CVSPGMR linear header file */
#include "sundials_dense.h"        /* use generic DENSE linear solver */

#ifdef _OPENMP
#include "nvector_openmp.h"	/* contains the definition of type N_Vector for openmp  */
#else
#include "nvector_serial.h"	/* contains the definition of type N_Vector  for serial */
#endif
#include "sundials_math.h"      /* contains UnitRoundoff, RSqrt,
                                 * SQR functions  */
#include "cvode_dense.h"        /* CVDENSE header file */

/*
 * SUNDIAL Header Files
 */
/* CVODE header file */
#include "cvode.h"
/* CVSPGMR linear header file */
#include "cvode_spgmr.h"
/* Definition of type N_Vector */
#ifdef _CVODE_OMP
#include "nvector_openmp.h"
#else
#include "nvector_serial.h"
#endif
/* UnitRoundoff, RSqrt, SQR functions */
#include "sundials_math.h"
/* CVDENSE header file */
#include "cvode_dense.h"

#ifdef _NOAH_
#include "spa.h"
#endif

#include "pihm_const.h"
#include "pihm_input_struct.h"
#include "elem_struct.h"
#include "river_struct.h"
#include "pihm_struct.h"
#include "pihm_func.h"
#endif
