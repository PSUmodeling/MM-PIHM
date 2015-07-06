#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

#define MULTF		2.0
#define PSIMIN		-70.0
#define EPS		0.05
#define THRESH		0.0
#define GRAV		9.80665
#define PI		3.14159265
#define BADVAL		-999
#define MAXSTRING	1024

#define MAT_CTRL        0
#define APP_CTRL        1
#define MAC_CTRL        2
#define SAT_CTRL        3

#define RIGHT_SIDE      0
#define LEFT_SIDE       1

#define CP           1004.0
#define LVH2O              2.503e6
#define SIGMA           5.67e-8
#define RD           287.04
#define RV             461.5

#define NUM_TS          5
#define NUM_METEO_TS    7
#define NUM_PRINT       100

extern int      verbose_mode;
extern int      debug_mode;

/* Enumrate type for forcing time series */
enum meteo_forcing_type {PRCP_TS, SFCTMP_TS, RH_TS, SFCSPD_TS, SOLAR_TS, LONGWAVE_TS, PRES_TS};
enum forcing_type {RIV_TS, BC_TS, METEO_TS, LAI_TS, SS_TS};

enum pihm_print_type {GW_CTRL, SURF_CTRL, SNOW_CTRL, RIVSTG_CTRL,
    INFIL_CTRL, RECHARGE_CTRL, CMC_CTRL, UNSAT_CTRL, EC_CTRL, ETT_CTRL, EDIR_CTRL,
    RIVFLX0_CTRL, RIVFLX1_CTRL, RIVFLX2_CTRL, RIVFLX3_CTRL, RIVFLX4_CTRL,
    RIVFLX5_CTRL, RIVFLX6_CTRL, RIVFLX7_CTRL, RIVFLX8_CTRL, RIVFLX9_CTRL,
    RIVFLX10_CTRL, SUBFLX_CTRL};
#endif
