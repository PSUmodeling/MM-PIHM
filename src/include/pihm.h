#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>

/* SUNDIAL Header Files */
#include "sundials_types.h"     /* realtype, integertype, booleantype
                                 * defination */
#include "cvode.h"              /* CVODE header file */
#include "cvode_spgmr.h"        /* CVSPGMR linear header file */
#include "sundials_smalldense.h"        /* use generic DENSE linear solver
                                         * for "small" */
#include "nvector_serial.h"     /* contains the definition of type
                                 * N_Vector */
#include "sundials_math.h"      /* contains UnitRoundoff, RSqrt,
                                 * SQR functions  */
#include "cvode_dense.h"        /* CVDENSE header file */

#ifdef _ENKF_
#include "mpi.h"
#endif
#ifdef _NOAH_
#include "spa.h"
#endif

#include "pihm_const.h"
#ifdef _ENKF_
#include "enkf.h"
#endif
#include "pihm_input_struct.h"
#include "elem_struct.h"
#include "river_struct.h"
#include "pihm_struct.h"
#include "pihm_func.h"
#endif
