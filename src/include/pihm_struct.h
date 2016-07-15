#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

#ifdef _BGC_
/* a structure to hold information on the annual co2 concentration */
typedef struct
{
    int             varco2;     /* (flag) 0=const 1=use file 2=const,file for Ndep */
    double          co2ppm;     /* (ppm)  constant CO2 concentration */
    double         *co2ppm_array;       /* (ppm)  annual CO2 concentration array */
    int            *co2year_array;      /* (year) year corresponding to the concentration value in co2ppm_arry */
    int             co2vals;    /* (num)  The number of CO2 concentration values in the co2ppm_array */
} co2control_struct;

/* a structure to hold annual nitrogen deposition data */
typedef struct
{
    int             varndep;    /* (flag) 0=const 1=use file  */
    double         *ndep_array; /* (kgN m-2 yr-1)  annual ndep array */
    int            *ndepyear_array;     /* (year) year corresponding to the ndep value in ndep_array */
    int             ndepvals;   /* (num)  The number of ndep values in the ndep_array */
    double          ndep;       /* (kgN/m2/yr) wet+dry atmospheric deposition of N */
    double          nfix;       /* (kgN/m2/yr) symbiotic+asymbiotic fixation of N */
} ndepcontrol_struct;

/* meteorological variable arrays */
#endif

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
#ifdef _BGC_
    double          simstarttime;       /* start time of simulation */
    double          simendtime; /* end time of simulation */
    int             bgc_spinup;     /* (flag) 1=spinup run, 0=normal run */
    int             maxspinyears;       /* maximum number of years for spinup run */
    int             dodaily;    /* flag for daily output */
    int             domonavg;   /* flag for monthly average of daily outputs */
    int             doannavg;   /* flag for annual average of daily outputs */
    int             doannual;   /* flag for annual output */
    int             ndayout;    /* number of daily outputs */
    int             nannout;    /* number of annual outputs */
    int            *daycodes;   /* array of indices for daily outputs */
    int            *anncodes;   /* array of indices for annual outputs */
    int             read_restart;       /* flag to read restart file */
    int             write_restart;      /* flag to write restart file */
    int             keep_metyr; /* (flag) 1=retain restart metyr, 0=reset metyr */
    int             onscreen;   /* (flag) 1=show progress on-screen 0=don't */
    int             spinupstartyear;    /* first met year for spinup */
    int             spinupendyear;      /* last met year for spinup */
    int             spinupstart;        /* start time of spinup */
    int             spinupend;  /* end time of spinup */
    cstate_struct   cs;
    nstate_struct   ns;
    cinit_struct    cinit;
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
#ifdef _CYCLES_
    agtbl_struct    agtbl;
    croptbl_struct  croptbl;
    mgmttbl_struct  mgmttbl;
#endif
#ifdef _BGC_
    co2control_struct co2;      /* CO2 concentration information */
    ndepcontrol_struct ndepctrl;        /* Nitrogen deposition control structure */
    epctbl_struct   epctbl;
#endif
    forc_struct     forc;

    elem_struct    *elem;       /* Store Element Information */
    river_struct   *riv;

    calib_struct    cal;
    ctrl_struct     ctrl;
    prtctrl_struct  prtctrl[NUM_PRINT];
}              *pihm_struct;

#endif
