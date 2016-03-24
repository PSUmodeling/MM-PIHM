#ifndef _ENKF_STRUCT_HEADER_
#define _ENKF_STRUCT_HEADER_
typedef struct var_struct
{
    char            name[MAXSTRING];
    int             dim;
    double          min;
    double          max;
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
    int             type;
    double          x;
    double          y;
    double          rad;
    double          depth;
    int             dim;
    int             nlyr;
    int            *var_ind;
    /* The linear observation operator vector has a general form of
     *
     * xf = SUM_j { SUM_i  [w_i * (k_ij x_ij + b_ij)]}
     *
     * w_i, k_j, and b_j are initialized at the beginning and stored here for
     * easy access so they don't need to be calculated at every EnKF step.
     * w_i could be the area of each model grid
     * k_j is ususally 1.0, or could be the thickness of each soil layer
     * b_j could be the soil column depth if xf is water table depth. */
    double         *weight;
    double        **k;
    double        **b;
} obs_struct;

typedef struct ensmbr_struct
{
    double          param[MAXPARAM];
    double         *var[MAXVAR];
} ensmbr_struct;

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
    int             ascii;
    int             numele;
    int             numriv;
    double          weight;
    ensmbr_struct  *member;
    int             update_var;
    int             update_param;

    //double         *total_water_volume;

    param_struct    param[MAXPARAM];
    var_struct      var[MAXVAR];
    obs_struct     *obs;
}              *enkf_struct;
#endif
