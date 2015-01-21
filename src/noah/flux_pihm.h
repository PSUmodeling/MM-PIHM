#ifndef FLUX_PIHM_HEADER
#define FLUX_PIHM_HEADER

#include "../pihm.h"
#include "noah.h"

void            LSM_initialize (char *, Model_Data, Control_Data *,
   LSM_STRUCT);
void            LSM_read (char *, LSM_STRUCT);
void            LSM_initialize_output (char *, Model_Data, LSM_STRUCT,
   char *);
void            LSM_PrintInit (Model_Data, LSM_STRUCT, char *);
void            LSM_FreeData (Model_Data, LSM_STRUCT);

void            PIHM2Noah (realtype, realtype, Model_Data, LSM_STRUCT);
void            Noah2PIHM (Model_Data, LSM_STRUCT);

int             FindLayer (LSM_STRUCT, double);
double          mod (double, double);

double topo_radiation (double Sdir, double Sdif, double zenith, double azimuth180, double slope, double aspect, double *h_phi, double svf);
#endif
