#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

/* 
 * Topographic parameters
 */
typedef struct topo_struct
{
    double          area;       /* area of element */
    double          x;          /* x of centroid */
    double          y;          /* y of centroid */
    double          zmin;       /* z_min of centroid */
    double          zmax;       /* z_max of centroid */
    double          zbed;
    double          edge[3];    /* edge i is from node i to node i+1 */
    double          surfx[3];
    double          surfy[3];
    double          node_zmax;
#ifdef _NOAH_
    double          slope;
    double          aspect;     /* ys: surface aspect of grid */
    double          svf;        /* ys: sky view factor */
    double          h_phi[36];  /* unobstrcted angle */
#endif
#ifdef _RT_
    double          areasub[3];
#endif
} topo_struct;

/*
 * Soil parameters
 */
typedef struct soil_struct
{
    double          depth;
    double          ksath;      /* horizontal geologic saturated
                                 * hydraulic conductivity */
    double          ksatv;      /* vertical geologic saturated
                                 * hydraulic conductivity */
    double          kinfv;      /* vertical surface saturated hydraulic
                                 * conductivity */
    double          dinf;       /* depth from ground surface accross
                                 * which head is calculated during
                                 *  infiltration */
    double          alpha;      /* alpha from van Genuchten eqn */
    double          beta;
    double          porosity;
    double          smcmax;
    double          smcmin;
    double          smcwlt;     /* wilting point (volumetric) */
    double          smcref;     /* soil moisture threshold where transpiration begins to stress (volumetric) */
    double          dmac;       /* macropore Depth */
    double          kmach;      /* macropore horizontal saturated
                                 * hydraulic conductivity */
    double          kmacv;      /* macropore vertical saturated
                                 * hydraulic conductivity */
    double          areafv;     /* macropore area fraction on a
                                 * vertical cross-section */
    double          areafh;     /* macropore area fraction on a
                                 * horizontal cross-section */
    int             macropore;  /* 1: macropore; 0: regular soil */
#ifdef _NOAH_
    double          csoil;      /* soil heat capacity (j m-3 k-1) */
    double          quartz;     /* soil quartz content */
    double          smcdry;     /* dry soil moisture threshold where direct evap frm top layer ends (volumetric) */
#endif
} soil_struct;

/*
 * Land cover ecophysiological paraemters
 */
typedef struct lc_struct
{
    int             type;
    double          shdfac;     /* areal fractional coverage of green vegetation (fraction= 0.0-1.0) */
    double          shdmin;     /* minimum areal fractional coverage of green vegetation (fraction= 0.0-1.0) <= shdfac */
    double          shdmax;     /* maximum areal fractional coverage of green vegetation (fraction= 0.0-1.0) >= shdfac */
    double          rzd;
    double          rsmin;      /* minimum canopy resistance */
    double          rgl;        /* reference incoming solar flux for
                                 * photosynthetically active canopy */
    double          laimin;
    double          laimax;     /* maxm. LAI accross all seasons for a
                                 * vegetation type */
    double          snup;       /* threshold snow depth (in water equivalent m) that implies 100 percent snow cover */
    double          hs;         /* parameter used in vapor pressure deficit function */
    double          topt;       /* optimum transpiration air temperature */
    double          cfactr;     /* parameter used in the canopy inteception calculation */
    double          rsmax;      /* max. stomatal resistance */
    double          emissmax;   /* maximum emissivity */
    double          emissmin;   /* minimum emmisivity */
    double          albedomax;  /* maximum background albedo */
    double          albedomin;  /* manimum background albedo */
    double          z0max;      /* maximum roughness length (m) */
    double          z0min;      /* minimum roughness length (m) */
    double          rough;      /* surface roughness of an element */
    double          cmcfactr;
    int             bare;
#ifdef _NOAH_
    int             isurban;
    int             nroot;      /* number of root layers, a function of veg type, determined in subroutine redprm. */
    double          ptu;        /* photo thermal unit (plant phenology for annuals/crops) (not yet used, but passed to redprm for future use in veg parms) */
    double          rtdis[MAXLYR];
#endif
} lc_struct;

/*
 * Land cover physical states
 */
typedef struct ps_struct
{
    int             macpore_status;
    double          rc;         /* lcp: canopy resistance (s m-1) */
    double          pc;         /* plant coefficient (unitless fraction, 0-1) where pc*etp = actual transp */
    double          xlai;       /* lcp: leaf area index (dimensionless) */
    double          rcs;        /* incoming solar rc factor (dimensionless) */
    double          rct;        /* air temperature rc factor (dimensionless) */
    double          rcq;        /* atmos vapor pressure deficit rc factor (dimensionless) */
    double          rcsoil;     /* soil moisture rc factor (dimensionless) */
    double          alb;        /* backround snow-free surface albedo (fraction), for julian day of year (usually from temporal interpolation of monthly */
    double          albedo;
    double          zlvl;       /* lcp: height (m) above ground of atmospheric forcing
                                 * variables */
    double          zlvl_wind;  /* lcp: height (m) above ground of wind observations */
    double          sfcspd;
    double          rh;
    double          sfcprs;
    double          surfavail;
#ifdef _NOAH_
    double          snoalb;     /* upper bound on maximum albedo over deep snow (e.g. from robinson and kukla, 1985, j. clim. & appl. meteor.) */
    int             nsoil;      /* number of soil layers (at least 2, and not
                                 * greater than parameter nsold set below) */

    double          sldpth[MAXLYR];     /* the thickness of each soil layer (m) */
    double          frzk;       /* frozen ground parameter */
    double          frzx;
    double          czil;       /* calculate roughness length of heat */
    double          emissi;     /* lcp: surface emissivity (between 0 and 1) */
    double          ch;         /* lcp: surface exchange coefficient for heat and moisture (m s-1); */
    double          rch;
    /*
     * note: ch is technically a conductance since it has been multiplied by wind speed. 
     */
    double          cm;         /* lcp: surface exchange coefficient for momentum (m s-1); */
    /*
     * note: cm is technically a conductance since it has been multiplied by wind speed. 
     */
    double          z0;         /* lcp: time varying roughness length (m) as function of snow depth */
    double          fcr;        /* YS: reduction of infiltration caused
                                 * by frozen ground */
    int             nmacd;
    double          salp;       /* shape parameter of distribution function of snow cover */
    double          fxexp;      /* soil evaporation exponent used in devap */
    double          sbeta;      /* parameter used to calculate vegetation effect on soil heat */
    double          lvcoef;
    double          snotime1;   /* age of the snow on the ground */
    double          ribb;       /* bulk richardson number used to limit the dew/frost */
    double          beta;       /* ratio of actual/potential evap (dimensionless) */
    double          sncovr;     /* fractional snow cover (unitless fraction, 0-1) */
    double          q1;         /* effective mixing ratio at surface (kg kg-1), used for diagnosing the mixing ratio at 2 meter for coupled model */
    double          q2;         /* mixing ratio at height zlvl above ground (kg kg-1) */
    double          cosz;       /* solar zenith angle (not used for now) */
    double          ffrozp;     /* fraction of frozen precipitation */
    double          z0brd;      /* background fixed roughness length (m) */
    double          embrd;      /* background surface emissivity (between 0 and 1) */
    double          q2sat;
    double          dqsdt2;
    int             nwtbl;
    double          sndens;
    double          snowh;      /* actual snow depth (m) */
    double          sncond;
    double          rr;         /* (-) */
    double          eta_kinematic;      /* atctual latent heat flux in kg m-2 s-1 */
    double          zbot;       /* depth (m) of lower boundary soil temperature */
    double          tbot;       /* bottom soil temperature (local yearly-mean sfc air temperature) */
    double          gwet;
    double          satdpth[MAXLYR];
#endif
} ps_struct;

typedef struct elemforc_struct
{
    int             bc_type[3];
    int             lai_type;

    /* Forcing values */
    double         *bc[3];
    double         *meteo[NUM_METEO_VAR];
    double         *lai;
    double         *z0;
    double         *source;
    double         *meltf;
    double         *riverbc;
#ifdef _NOAH_
    double         *rad[2];
#endif
} elemforc_struct;

#ifdef _DAILY_
typedef struct daily_struct
{
    double          sfctmp;
    double          tday;
    double          tnight;
    double          tmax;
    double          tmin;

    double          solar;
    double          solar_total;

    double          sfcprs;

    double          surf;
    double          unsat;
    double          gw;

    double          fluxsurf[4];
    double          fluxsub[4];
    double          infil;
    double          rechg;

    double          dayl;
    double          prev_dayl;

    int             counter;
    int             daylight_counter;

#ifdef _NOAH_
    double          stc;
    double          sh2o;
    double          q2d;
    double          albedo;
    double          ch;
#endif
} daily_struct;
#endif

typedef struct ws_struct
{
    double          stage;
    double          surf;
    double          gw;
    double          unsat;
    double          sneqv;      /* liquid water-equivalent snow depth (m). note: snow density = sneqv/snowh */
    double          cmcmax;     /* maximum canopy water capacity */
    double          cmc;      /* Interception storage */
#ifdef _NOAH_
    double          smc[MAXLYR];        /* total soil moisture content (volumetric fraction) */
    double          sh2o[MAXLYR];       /* unfrozen soil moisture content (volumetric fraction). note: frozen soil moisture = smc - sh2o */
    double          soilw;      /* available soil moisture in root zone (unitless fraction between smcwlt and smcmax) */
    double          soilm;      /* total soil column moisture content (frozen+unfrozen) (m) */
#endif
} ws_struct;

typedef struct wf_struct
{
    double          fluxsurf[3];        /* Overland Flux */
    double          fluxsub[3]; /* Subsurface Flux */
    double          prcp;       /* Precep. on each element */
    double          netprcp;    /* Net precep. on each elment */
    double          infil;      /* Variable infiltration rate */
    double          rechg;      /* Recharge rate to GW */
    double          drip;       /* through-fall of precip and/or dew in excess of canopy water-holding capacity (m/s) */
    double          edir;
    double          ett;
    double          ec;
    double          etp;
    double          eta;
    double          et[MAXLYR];
    double          edir_sfc;
    double          edir_unsat;
    double          edir_gw;
    double          ett_unsat;
    double          ett_gw;
    double          fluxriv[11];
#ifdef _NOAH_
    double          runoff1;    /* surface runoff (m s-1), not infiltrating the surface */
    double          runoff2;    /* subsurface runoff (m s-1), drainage out bottom of last soil layer (baseflow) */
    double          runoff2_lyr[MAXLYR];
    double          runoff3;    /* numerical trunctation in excess of porosity (smcmax) for a given soil layer at the end of a time step (m s-1). note: the above runoff2 is actually the sum of runoff2 and runoff3 */
    double          pcpdrp;     /* combined prcp1 and drip (from cmc) that goes into the soil (m s-1) */
    double          prcprain;   /* liquid-precipitation rate (kg m-2 s-1) (not used) */
    double          dew;        /* dewfall (or frostfall for t<273.15) (m) */
    double          snomlt;     /* snow melt (m/s) (water equivalent) */
    double          esnow;      /* sublimation from (or deposition to if <0) snowpack (ms-1) */
    double          etns;
#endif
} wf_struct;

typedef struct es_struct
{
    double          stc[MAXLYR];        /* soil temp (k) */
    double          t1;         /* ground/canopy/snowpack) effective skin temperature (k) */
    double          th2;        /* air potential temperature (k) at height zlvl above ground */
    double          sfctmp;
} es_struct;

typedef struct ef_struct
{
    double          solnet;     /* net downward solar radiation ((w m-2; positive) */
    double          etp;        /* potential evaporation (w m-2) */
    double          epsca;      /* */
    double          ssoil;      /* soil heat flux (w m-2: negative if downward from surface) */
    double          eta;        /* actual latent heat flux (w m-2: negative, if up from surface) */
    double          sheat;      /* sensible heat flux (w m-2: negative, if upward from surface) */
    double          fdown;      /* radiation forcing at the surface (w m-2) = soldn*(1-alb)+lwdn */
    double          lwdn;
    double          ec;         /* canopy water evaporation (w m-2) */
    double          edir;       /* direct soil evaporation (w m-2) */
    double          et[MAXLYR]; /* plant transpiration from a particular root (soil) layer (w m-2) */
    double          ett;        /* total plant transpiration (w m-2) */
    double          esnow;      /* sublimation from (or deposition to if <0) snowpack (w m-2) */
    double          soldn;      /* solar downward radiation (w m-2; positive, not net solar) */
    double          longwave;
    double          flx1;       /* precip-snow sfc (w m-2) */
    double          flx2;       /* freezing rain latent heat flux (w m-2) */
    double          flx3;       /* phase-change heat flux from snowmelt (w m-2) */
    double          solardirect;        /* direct component of downward solar radiation (w m-2) (not used) */
} ef_struct;

typedef struct elemic_struct
{
    double          surf;
    double          unsat;
    double          gw;
    double          sneqv;
    double          intcp;
#ifdef _NOAH_
    double          t1;
    double          snowh;
    double          stc[MAXLYR];
    double          smc[MAXLYR];
    double          sh2o[MAXLYR];
#endif
} elemic_struct;

typedef struct elem_struct
{
    int             node[3];    /* Counterclock-wise */
    int             nabr[3];    /* neighbor i shares edge i
                                 * (0: on boundary) */
    int             ind;

    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
    elemforc_struct forc;
    elemic_struct   ic;

#ifdef _DAILY_
    daily_struct    daily;
#endif
    ps_struct       ps;

    ws_struct       ws;
    ws_struct       ws0;
    wf_struct       wf;
#ifdef _NOAH_
    wf_struct       avgwf;

    es_struct       es;
    ef_struct       ef;
#endif

} elem_struct;

typedef struct shp_struct
{
    double          depth;
    int             intrpl_ord;
    double          coeff;
    double          length;
} shp_struct;

typedef struct matl_struct
{
    double          rough;
    double          cwr;
    double          ksath;
    double          ksatv;
    double          bedthick;
    double          porosity;
} matl_struct;

typedef struct riveric_struct
{
    double          stage;
    double          gw;
} riveric_struct;

typedef struct river_struct
{
    topo_struct     topo;
    shp_struct      shp;
    matl_struct     matl;
    elemforc_struct forc;
    ws_struct       ws;
    ws_struct       ws0;
    wf_struct       wf;
    riveric_struct  ic;
#ifdef _DAILY_
    daily_struct    daily;
#endif
    int             leftele;    /* Left neighboring element */
    int             rightele;   /* Right neighboring element */
    int             fromnode;   /* Upstream Node no. */
    int             tonode;     /* Dnstream Node no. */
    int             down;       /* down stream segment */
} river_struct;

typedef struct calib_struct
{
    double          ksath;
    double          ksatv;
    double          kinfv;
    double          kmach;
    double          kmacv;
    double          dinf;
    double          rzd;
    double          dmac;
    double          porosity;
    double          alpha;
    double          beta;
    double          areafv;
    double          areafh;
    double          vegfrac;
    double          albedo;
    double          rough;
    double          prcp;
    double          sfctmp;

    double          ec;
    double          ett;
    double          edir;

    double          rivrough;
    double          rivksath;
    double          rivksatv;
    double          rivbedthick;
    double          rivdepth;
    double          rivshpcoeff;

#ifdef _NOAH_
    double          thetaref;   /* ys */
    double          thetaw;     /* ys */
    double          rsmin;      /* ys */
    double          drip;       /* ys */
    double          intcp;      /* ys */
    double          czil;       /* ys */
    double          fxexp;      /* ys */
    double          cfactr;     /* ys */
    double          rgl;        /* ys */
    double          hs;         /* ys */
#endif
#ifdef _RT_
    double          pco2;
    double          keq;
    double          ssa;
    double          site_den;
    double          prep_conc;
#endif
} calib_struct;

typedef struct ctrl_struct
{
    int             ascii;      /* ys: add ascii output model
                                 * (default is binary */
    int             write_ic;   /* ys: runs model as spinup. model output at
                                 * the last step will be saved in .init */
    int             solver;     /* solver type */
    int             nstep;      /* number of external time steps (when
                                 * results can be printed) for the
                                 * whole simulation */
    int             nprint;     /* ys: number of variables for output */

    /* time interval to output average values of variables
     * variables will not be printed if time intervals are set to 0 */
    int             prtvrbl[NUM_PRINT];

    int             init_type;  /* initialization mode */

    int             unsat_mode; /* unsat mode */
    int             surf_mode;  /* surface overland flow mode */
    int             riv_mode;   /* river routing mode */


    double          abstol;     /* absolute tolerance */
    double          reltol;     /* relative tolerance */
    double          initstep;   /* initial step size */
    double          maxstep;    /* maximum step size */

    int             etstep;     /* step for et from interception */

    int             starttime;  /* start time of simulation */
    int             endtime;    /* end time of simulation */

    int             stepsize;

    int            *tout;
#ifdef _NOAH_
    int             nsoil;
    double          sldpth[MAXLYR];
    int             rad_mode;   /* radiation mode; 1: topographic, 0: uniform */
#endif
} ctrl_struct;

typedef struct prtctrl_struct
{
    char            name[MAXSTRING];
    int             intvl;
    int             nvrbl;
    double        **vrbl;
    double         *buffer;
} prtctrl_struct;

/* Model_data definition */
typedef struct pihm_struct
{
    int             numele;     /* number of elements */
    int             numriv;     /* number of rivere segments */

    double          longitude;
    double          latitude;
    double          elevation;

    /*
     * Input
     */
    /* File names */
    filename_struct filename;

    /* Input look-up talbes */
    meshtbl_struct  meshtbl;
    atttbl_struct   atttbl;
    soiltbl_struct  soiltbl;    /* Store Soil Information */
    geoltbl_struct  geoltbl;    /* Store Soil Information */
    lctbl_struct    lctbl;      /* Store Land Cover Information */
    rivtbl_struct   rivtbl;     /* Store River Segment Information */
    shptbl_struct   shptbl;     /* Store River Shape Information */
    matltbl_struct  matltbl;    /* Store River Bank Material Information */
#ifdef _NOAH_
    noahtbl_struct  noahtbl;
#endif
    forc_struct     forc;

    elem_struct    *elem;       /* Store Element Information */
    river_struct   *riv;

    calib_struct    cal;
    ctrl_struct     ctrl;
    prtctrl_struct  prtctrl[NUM_PRINT];
}              *pihm_struct;

#endif
