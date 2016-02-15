#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

typedef struct filename_struct
{
    char            riv[MAXSTRING];
    char            mesh[MAXSTRING];
    char            att[MAXSTRING];
    char            soil[MAXSTRING];
    char            geol[MAXSTRING];
    char            lc[MAXSTRING];
    char            forc[MAXSTRING];
    char            lai[MAXSTRING];
    char            ibc[MAXSTRING];
    char            para[MAXSTRING];
    char            calib[MAXSTRING];
    char            init[MAXSTRING];
#ifdef _NOAH_
    char            lsm[MAXSTRING];
    char            rad[MAXSTRING];
    char            lsminit[MAXSTRING];
#endif
} filename_struct;

typedef struct rivtbl_struct
{
    int             number;
    int            *fromnode;   /* Upstream Node no. */
    int            *tonode;     /* Dnstream Node no. */
    int            *down;       /* down stream segment */
    int            *leftele;    /* Left neighboring element */
    int            *rightele;   /* Right neighboring element */
    int            *shp;        /* shape type    */
    int            *matl;       /* material type */
    int            *bc;         /* BC type */
    int            *rsvr;
} rivtbl_struct;

typedef struct shptbl_struct
{
    int             number;
    double         *depth;      /* depth */
    int            *intrpl_ord; /* Interpolation order for river shape:
                                 * 1: rectangle,
                                 * 2: triangle,
                                 * 3: quadratic,
                                 * 4: cubic */
    double         *coeff;      /* Coefficient c in
                                 * D = c * pow(B / 2, interpOrd) */
} shptbl_struct;

typedef struct matltbl_struct
{
    int             number;     /* Number of River Bank/Bed Material */
    double         *rough;
    double         *cwr;        /* Weir Discharge Coefficient */
    double         *ksath;      /* Conductivity of river banks */
    double         *ksatv;      /* Conductivity of river bed */
    double         *bedthick;   /* thickeness of conductive river bed */
} matltbl_struct;

typedef struct meshtbl_struct
{
    int             numele;
    int             numnode;
    int           **node;
    int           **nabr;
    double         *x;
    double         *y;
    double         *zmin;
    double         *zmax;
} meshtbl_struct;

typedef struct atttbl_struct
{
    int            *soil;       /* soil type */
    int            *geol;
    int            *lc;         /* Land Cover type  */
    int           **bc;         /* Boundary condition type.
                                 * 0: Natural BC (no flow);
                                 * 1: Dirichlet BC;
                                 * 2:Neumann BC */
    int            *meteo;      /* precipitation (forcing) type */
    int            *lai;        /* LAI forcing type (0: use climatological
                                 * values; else: use user provided time
                                 *                                  * series */
    int            *source;     /* source (well) type */
    int            *macropore;
} atttbl_struct;

typedef struct soiltbl_struct
{
    int             number;     /* index */
    int            *mukey;
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *kinfv;
    double         *ksatv;      /* vertical saturated soil
                                 * conductivity */
    double         *ksath;
    double         *smcmax;     /* soil porosity */
    double         *smcmin;     /* soil moisture residual */
    double         *qtz;        /* ys: quartz content */
    double         *alpha;      /* soil curve parameter 1 */
    double         *beta;       /* soil curve parameter 2 */

    double         *areafh;     /* macroporous area fraction on
                                 * horizontal section */
    double         *areafv;
    double         *kmacv;      /* macroporous saturated vertical
                                 * conductivity */
    double         *kmach;
    double         *dmac;
    double         *dinf;       /* depth from ground surface accross which
                                 * head is calculated during infiltration */
    double         *smcref;
    double         *smcwlt;
} soiltbl_struct;

typedef struct geoltbl_struct
{
    int             number;     /* index */
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *ksath;      /* horizontal saturated geology
                                 * conductivity */
    double         *ksatv;      /* vertical saturated geology
                                 * conductivity */
    double         *smcmax;     /* geology porosity */
    double         *smcmin;    /* residual porosity */
    double         *alpha;      /* van genuchten parameter */
    double         *beta;       /* van genuchten parameter */
} geoltbl_struct;

typedef struct lctbl_struct
{
    int             number;     /* index */

    double         *laimax;     /* max lai */
    double         *laimin;     /* ys: min lai */
    double         *vegfrac;    /* canopy fracn */
    double         *albedomin;  /* ys: minimum albedo */
    double         *albedomax;  /* ys: maximum albedo */
    double         *emissmin;   /* ys: minimum emissivity */
    double         *emissmax;   /* ys: maximum emissivity */
    double         *z0min;      /* ys: minimum roughness length */
    double         *z0max;      /* ys: maximum roughness length */
    double         *hs;         /* ys: vapor pressure deficit stress
                                 * parameter */
    double         *snup;       /* ys */
    double         *rgl;        /* visible solar flux used in radiation
                                 * stress */
    double         *rsmin;      /* minimum stomatal resistance */
    double         *rough;      /* surface roughness factor  */
    double         *rzd;        /* rootzone depth */

    double          rsmax;      /* YS */
    int             bare;       /* YS */
    int             natural;
    double          cfactr;     /* YS */
    double          topt;       /* YS */
} lctbl_struct;

typedef struct tsdata_struct
{
    int             length;     /* length of time series */
    int            *ftime;
    double        **data;       /* 2D time series data */
    double         *value;
    double          zlvl_wind;
} tsdata_struct;

typedef struct forc_struct
{
    /* Forcing series */
    int             nbc;
    tsdata_struct  *bc;

    int             nmeteo;
    tsdata_struct  *meteo;

    int             nlai;
    tsdata_struct  *lai;

    int             nz0;
    tsdata_struct  *z0;

    int             nsource;
    tsdata_struct  *source;

    int             nmeltf;
    tsdata_struct  *meltf;

    int             nriverbc;
    tsdata_struct  *riverbc;
#ifdef _NOAH_
    int             nrad;
    tsdata_struct  *rad;
#endif
} forc_struct;

#ifdef _NOAH_
typedef struct noahtbl_struct
{
    double          sbeta;
    double          fxexp;
    double          csoil;
    double          salp;
    double          refdk;
    double          refkdt;
    double          frzk;
    double          zbot;
    double          tbot;
    double          smlow;
    double          smhigh;
    double          czil;
    double          lvcoef;
} noahtbl_struct;
#endif
#endif
