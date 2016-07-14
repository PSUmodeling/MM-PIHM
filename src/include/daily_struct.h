#ifndef _DAILY_STRUCT_
#define _DAILY_STRUCT_

typedef struct daily_wstate_struct
{
    double          avg_surf;
    double          avg_unsat;
    double          avg_gw;
    double          avg_sh2o[MAXLYR];
} daily_wstate_struct;

typedef struct daily_wflux_struct
{
    double          avg_surf[NUM_EDGE];
    double          avg_subsurf[NUM_EDGE];
    double          avg_et[MAXLYR];
    double          avg_smflxv[MAXLYR];
} daily_wflux_struct;

typedef struct daily_pstate_struct
{
    double          dayl;
    double          prev_dayl;
    double          avg_q2d;
    double          avg_sfcprs;
    double          avg_ch;
    double          avg_albedo;
    double          avg_sfcspd;
} daily_pstate_struct;

typedef struct daily_estate_struct
{
    double          tmax;
    double          tmin;
    double          avg_sfctmp;
    double          tday;
    double          tnight;
    double          avg_stc[MAXLYR];
} daily_estate_struct;

typedef struct daily_eflux_struct
{
    double          avg_soldn;
    double          par;
    double          solar_total;
} daily_eflux_struct;

typedef struct daily_struct
{
    int             counter;
    int             daylight_counter;
    daily_wstate_struct ws;
    daily_wflux_struct  wf;
    daily_pstate_struct ps;
    daily_estate_struct es;
    daily_eflux_struct  ef;
} daily_struct;

typedef struct river_daily_wstate_struct
{
    double          avg_stage;
    double          avg_gw;
} river_daily_wstate_struct;

typedef struct river_daily_wflux_struct
{
    double          avg_river[NUM_RIVFLX];
} river_daily_wflux_struct;

typedef struct river_daily_struct
{
    int             counter;

    river_daily_wstate_struct   ws;
    river_daily_wflux_struct    wf;
} river_daily_struct;
#endif
