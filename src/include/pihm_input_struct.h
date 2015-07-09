#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

typedef struct riv_att_tbl_struct
{
    int             number;
    int            *fromnode;   /* Upstream Node no. */
    int            *tonode;     /* Dnstream Node no. */
    int            *down;       /* down stream segment */
    int            *leftele;    /* Left neighboring element */
    int            *rightele;   /* Right neighboring element */
    int            *shp;        /* shape type    */
    int            *matl;       /* material type */
    int            *ic;         /* IC type */
    int            *bc;         /* BC type */
    int            *rsvr;
} riv_att_tbl_struct;

typedef struct riv_shp_tbl_struct
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
} riv_shp_tbl_struct;

typedef struct riv_matl_tbl_struct
{
    int             number;     /* Number of River Bank/Bed Material */
    double         *rough;
    double         *cwr;        /* Weir Discharge Coefficient */
    double         *ksath;      /* Conductivity of river banks */
    double         *ksatv;      /* Conductivity of river bed */
    double         *bedthick;   /* thickeness of conductive river bed */
} riv_matl_tbl_struct;

typedef struct riv_ic_tbl_struct
{
    int             number;
    double         *stage;
} riv_ic_tbl_struct;

typedef struct mesh_tbl_struct
{
    int             numele;
    int             numnode;
    int           **node;
    int           **nabr;
    double         *x;
    double         *y;
    double         *zmin;
    double         *zmax;
} mesh_tbl_struct;

typedef struct attrib_tbl_struct
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
} attrib_tbl_struct;

typedef struct soil_tbl_struct
{
    int             number;     /* index */
    double         *ksatv;      /* vertical saturated soil
                                 * conductivity */
    double         *thetas;     /* soil porosity */
    double         *thetar;     /* soil moisture residual */
    double         *qtz;        /* ys: quartz content */
    double         *alpha;      /* soil curve parameter 1 */
    double         *beta;       /* soil curve parameter 2 */

    double         *areafh;     /* macroporous area fraction on
                                 * horizontal section */
    double         *kmacv;      /* macroporous saturated vertical
                                 * conductivity */
    double         *dinf;       /* depth from ground surface accross which
                                 * head is calculated during infiltration */
} soil_tbl_struct;

typedef struct geol_tbl_struct
{
    int             number;     /* index */
    double         *ksath;      /* horizontal saturated geology
                                 * conductivity */
    double         *ksatv;      /* vertical saturated geology
                                 * conductivity */
    double         *thetas;     /* geology porosity */
    double         *thetar;     /* residual porosity */
    double         *alpha;      /* van genuchten parameter */
    double         *beta;       /* van genuchten parameter */

    double         *areafv;     /* macroporous area fraction on vertical
                                 * section */
    double         *kmach;      /* macroporous saturated
                                 * horizontal conductivity */
    double         *dmac;
} geol_tbl_struct;

typedef struct lc_tbl_struct
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
} lc_tbl_struct;

typedef struct ts_struct
{
    int             length;     /* length of time series */
    int            *ftime;
    double        **data;       /* 2D time series data */
} ts_struct;

typedef struct forcing_ts_struct
{
    int             nts[NUM_TS];

    ts_struct      *ts[NUM_TS];
    double         *zlvl_wind;  /* wind measurement height */

    double         *bc;
    double         *meteo[NUM_METEO_TS];
    double         *lai;
    double         *z0;
    double         *source;
    double         *meltf;
    double         *riverbc;
} forcing_ts_struct;
#endif
