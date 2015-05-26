#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

#define multF		2.0
#define MINpsi		-70.0
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

#define C_air           1004.0
#define Lv              2.503e6
#define SIGMA           5.67e-8
#define R_dry           287.04
#define R_v             461.5

/* Enumrate type for forcing time series */
enum forcing_type {PRCP_TS, SFCTMP_TS, RH_TS, SFCSPD_TS, SOLAR_TS, LONGWAVE_TS, PRES_TS, LAI_TS, RL_TS, MF_TS, SS_TS};

#endif
