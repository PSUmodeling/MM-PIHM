#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

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

#ifdef _RT_
    double          areasub[3];
#endif
} topo_struct;

typedef struct soil_struct
{
    double          depth;
    double          ksath;      /* horizontal geologic saturated
                                 * hydraulic conductivity */
    double          ksatv;      /* vertical geologic saturated
                                 * hydraulic conductivity */
    double          kinfv;      /* vertical surface saturated hydraulic
                                 * conductivity */
    double          porosity;
    double          dinf;       /* depth from ground surface accross
                                 * which head is calculated during
                                 *  infiltration */
    double          alpha;      /* alpha from van Genuchten eqn */
    double          beta;
    double          thetas;
    double          thetar;
    double          thetaref;   /* YS: Soil field capacity */
    double          thetaw;     /* YS: Soil wilting point */
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
} soil_struct;

typedef struct lc_struct
{
    int             type;
    double          vegfrac;    /* areal vegetation fraction in a
                                 * triangular element */
    double          rzd;
    double          rsmin;      /* minimum canopy resistance */
    double          rgl;        /* reference incoming solar flux for
                                 * photosynthetically active canopy */
    double          hs;
    double          snup;
    double          laimin;
    double          laimax;     /* maxm. LAI accross all seasons for a
                                 * vegetation type */
    double          emissmin;
    double          emissmax;
    double          albedo;
    double          albedomin;
    double          albedomax;
    double          z0min;
    double          z0max;
    double          rough;      /* surface roughness of an element */
    double          intcp_factr;
    double          rsmax;
    int             bare;
    double          cfactr;
    double          topt;
} lc_struct;

typedef struct forc_struct
{
    double          zlvl_wind;
    int             bc_type[3];
    int             lai_type;

    /* Forcing values */
    double         *bc[3];
    double         *meteo[NUM_METEO_TS];
    double         *lai;
    double         *z0;
    double         *source;
    double         *meltf;
    double         *riverbc;
} forc_struct;

typedef struct daily_struct
{
    double          sfctmp;
    double          tday;
    double          tnight;
    double          tmax;
    double          tmin;

    double          solar;

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

typedef struct elem_struct
{
    int             node[3];    /* Counterclock-wise */
    int             nabr[3];    /* neighbor i shares edge i
                                 * (0: on boundary) */
    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
    forc_struct     forc;

    daily_struct    daily;

    double          surf0;      /* YS: Stores surface water level of
                                 * last time step */
    double          gw0;        /* YS: Stores groundwater level of last time
                                 * step */
    double          unsat0;     /* YS: Stores unsaturated storage of
                                 * last time step */
    double          surf;
    double          gw;
    double          unsat;
    double          fluxsurf[3];        /* Overland Flux */
    double          fluxsub[3]; /* Subsurface Flux */
    double          runoff;
    double          prcp;       /* Precep. on each element */
    double          netprcp;    /* Net precep. on each elment */
    double          infil;      /* Variable infiltration rate */
    double          rechg;      /* Recharge rate to GW */
    double          snow;       /* YS: Snow water equivalent on each element */
    double          intcp;      /* Interception storage */
    double          drip;       /* Through Fall */
    double          et[3];      /* Evapotranspiration (from canopy, ground,
                                 * transpiration) */

    double          edir[3];
    double          ett[3];

    int             macpore_status;

    double          albedo;
#ifdef _NOAH_
    double          sfcsat;     /* YS: Surface saturation */
    double          et_from_sat;        /* YS: Fraction of Transpiration that
                                         * extracts from the saturated zone */
    double          fcr;        /* YS: reduction of infiltration caused
                                 * by frozen ground */
    double          totalw;
    double          mbc;
#endif

#ifdef _RT_
    double          temp;       /* temperature   */
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

typedef struct river_struct
{
    topo_struct     topo;
    shp_struct      shp;
    matl_struct     matl;
    forc_struct     forc;

    daily_struct    daily;

    int             leftele;    /* Left neighboring element */
    int             rightele;   /* Right neighboring element */
    int             fromnode;   /* Upstream Node no. */
    int             tonode;     /* Dnstream Node no. */
    int             down;       /* down stream segment */

    double          stage0;
    double          stage;
    double          gw0;
    double          gw;
    double          fluxriv[11];
    double          totalw;
    double          mbc;
} river_struct;

/*
 * Initial state variable conditions on each element
 */
typedef struct ic_struct
{
    double         *intcp;      /* Interception storage (Note all these
                                 * variables have dimension of L */
    double         *snow;       /* Snow depth */
    double         *surf;       /* Overland flow depth */
    double         *unsat;      /* unsaturated zone depth */
    double         *gw;         /* saturated zone depth */
    double         *rivgw;
    double         *stage;
} ic_struct;

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

    double          et[3];

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

    /* Input */
    mesh_tbl_struct mesh_tbl;
    attrib_tbl_struct attrib_tbl;
    ic_struct       ic;         /* Store Element Initial Condtion */
    soil_tbl_struct soil_tbl;   /* Store Soil Information */
    geol_tbl_struct geol_tbl;   /* Store Soil Information */
    lc_tbl_struct   lc_tbl;     /* Store Land Cover Information */
    forcing_ts_struct forcing;
    riv_att_tbl_struct riv_att_tbl;     /* Store River Segment Information */
    riv_shp_tbl_struct riv_shp_tbl;     /* Store River Shape Information */
    riv_matl_tbl_struct riv_matl_tbl;   /* Store River Bank Material Information */
    riv_ic_tbl_struct riv_ic_tbl;
    //riv_rsvr_tbl_struct riv_rsvr_tbl;

    elem_struct    *elem;       /* Store Element Information */
    river_struct   *riv;

    calib_struct    cal;
    ctrl_struct     ctrl;
    prtctrl_struct  prtctrl[NUM_PRINT];
}              *pihm_struct;

#endif
