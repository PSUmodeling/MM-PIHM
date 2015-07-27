//#include "mpi.h"

#define MAXPARAM        100
#define MAXVAR          100

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
    double         *parameter;
    double        **variable;
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
    int             update_prmt;

    double         *total_water_volume;

    param_struct    param[MAXPARAM];
    var_struct      var[MAXVAR];
    obs_struct     *obs;
} *enkf_struct;

void EnKFRead(char *project, enkf_struct ens);

