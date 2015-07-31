#include "mpi.h"

#define MAXPARAM        100
#define MAXVAR          100
#define MAXINT          2147483647
#define CORRMAX         0.25

#define SUCCESS_TAG     2
#define NE_TAG          1

enum prmt_type {KSATH, KSATV, KINF, KMACH, KMACV, DINF, RZD, DMAC, POROSITY,
    ALPHA, BETA, AREAFV, AREAFH, VEGFRAC, ALBEDO, ROUGH, PRCP, SFCTMP,
    ET0, ET1, ET2, RIVROUGH, RIVKSATH, RIVKSATV, RIVBEDTHICK, RIVDEPTH,
    RIVSHPCOEFF, THETAREF, THETAW, RSMIN, DRIP, INTCP, CZIL, FXEXP, CFACTR,
    RGL, HS};

typedef struct var_struct
{
    char            name[MAXSTRING];
    int             dim;
} var_struct;

typedef struct param_struct
{
    int             perturb;
    int             update;
    double          perturb_min;
    double          perturb_max;
    double          init_std;
    double          min;
    double          max;
    int             type;
    char            name[MAXSTRING];
} param_struct;

typedef struct obs_struct
{
    char            name[MAXSTRING];
    char            fn[MAXSTRING];
    int             obs_type;
    int             nctrl;
    int            *var_ind;
    int            *grid_ind;
    double         *weight;
} obs_struct;

typedef struct ens_mbr_struct
{
    double         *param;
    double        **var;
} ens_mbr_struct;

typedef struct enkf_struct
{
    int             ne;
    int             nobs;
    int             interval;
    int             start_mode;
    int             end_time;
    int             cycle_start_time;
    int             cycle_end_time;
    int             mbr_start_mode;
    int             numele;
    int             numriv;
    double          weight;
    ens_mbr_struct *member;
    int             update_var;
    int             update_param;

    double         *total_water_volume;

    param_struct    param[MAXPARAM];
    var_struct      var[MAXVAR];
    obs_struct     *obs;
} *enkf_struct;

void EnKFRead(char *project, enkf_struct ens);
double Randn();
void MapVar (var_struct *var, int numele, int numriv);
void Perturb(char *project, enkf_struct ens, char *outputdir);
void MapVar (var_struct *var, int numele, int numriv);
void Calib2Mbr (calib_struct cal, double *param);
void WriteParamOutput (int rawtime, enkf_struct ens, int ind, char *outputdir);
void WriteCalFile (enkf_struct ens, char *project);
void JobHandout (int msg, int total_jobs);
void PrintEnKFStatus (int starttime, int endtime);
void JobHandIn (int total_jobs);
void WritePara (char *project, int start_mode, int start_time, int end_time);

