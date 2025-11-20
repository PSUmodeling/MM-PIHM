/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for a custom N_Vector implementation
 * (note that this template just implements the serial N_Vector).
 * -----------------------------------------------------------------*/

#ifndef _MY_NVECTOR_H
#define _MY_NVECTOR_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Custom implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _My_NVectorContent
{
  sunindextype length;     /* vector length       */
  sunbooleantype own_data; /* data ownership flag */
  sunrealtype* data;       /* data array          */
};

typedef struct _My_NVectorContent* My_NVectorContent;

/*
 * -----------------------------------------------------------------
 * Macros MY_NV_CONTENT, MY_NV_DATA, MY_NV_OWN_DATA,
 *        MY_NV_LENGTH, and MY_NV_Ith
 * -----------------------------------------------------------------
 */

#define MYNV_CONTENT(v) ((My_NVectorContent)(v->content))

#define MYNV_LENGTH(v) (MYNV_CONTENT(v)->length)

#define MYNV_OWN_DATA(v) (MYNV_CONTENT(v)->own_data)

#define MYNV_DATA(v) (MYNV_CONTENT(v)->data)

#define MYNV_Ith(v, i) (MYNV_DATA(v)[i])

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_serial
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
N_Vector N_VNewEmpty_Mine(sunindextype vec_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNew_Mine(sunindextype vec_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VMake_Mine(sunindextype vec_length, sunrealtype* v_data,
                      SUNContext sunctx);

SUNDIALS_EXPORT
sunindextype N_VGetLength_Mine(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_Mine(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_Mine(N_Vector v, FILE* outfile);

static inline N_Vector_ID N_VGetVectorID_Mine(N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_Mine(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_Mine(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_Mine(N_Vector v);

SUNDIALS_EXPORT
void N_VSpace_Mine(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
sunrealtype* N_VGetArrayPointer_Mine(N_Vector v);

SUNDIALS_EXPORT
void N_VSetArrayPointer_Mine(sunrealtype* v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT
void N_VLinearSum_Mine(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z);
SUNDIALS_EXPORT
void N_VConst_Mine(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_Mine(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_Mine(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_Mine(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_Mine(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_Mine(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_Mine(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VDotProd_Mine(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNorm_Mine(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_Mine(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_Mine(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VMin_Mine(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_Mine(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VL1Norm_Mine(N_Vector x);

SUNDIALS_EXPORT
void N_VCompare_Mine(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VInvTest_Mine(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMask_Mine(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotient_Mine(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_Mine(int nvec, sunrealtype* c, N_Vector* V,
                                     N_Vector z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_Mine(int nvec, sunrealtype* a, N_Vector x,
                                 N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_Mine(int nvec, N_Vector x, N_Vector* Y,
                                sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_Mine(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_Mine(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_Mine(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_Mine(int nvecs, N_Vector* X, N_Vector* W,
                                       sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_Mine(int nvecs, N_Vector* X, N_Vector* W,
                                           N_Vector id, sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMultiVectorArray_Mine(int nvec, int nsum, sunrealtype* a,
                                            N_Vector* X, N_Vector** Y,
                                            N_Vector** Z);

SUNDIALS_EXPORT
SUNErrCode N_VLinearCombinationVectorArray_Mine(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_Mine(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_Mine(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT
SUNErrCode N_VBufSize_Mine(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_Mine(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_Mine(N_Vector x, void* buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Mine(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Mine(N_Vector v,
                                                      sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
