#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/* Attribute */
typedef struct river_attrib_struct
{
    int             riverbc_type;
} attrib_struct;

/* 
 * Topographic parameters
 */
typedef struct river_topo_struct
{
    double          area;       /* area of element */
    double          x;          /* x of centroid */
    double          y;          /* y of centroid */
    double          zmin;       /* z_min of centroid */
    double          zmax;       /* z_max of centroid */
    double          zbed;
    double          node_zmax;
} river_topo_struct;

typedef struct river_wstate_struct
{
    double          stage;
    double          gw;
} river_wstate_struct;

typedef struct river_wflux_struct
{
    double          river[11];
} river_wflux_struct;

typedef struct river_stor_struct
{
    double         *stage;
    double         *gw;
    double         *riverflx[11];
    int            *flag;
} river_stor_struct;

#ifdef _DAILY_
typedef struct river_daily_struct
{
    int             counter;
    int             daylight_counter;

    river_wstate_struct ws;
    river_wflux_struct  wf;
} river_daily_struct;
#endif

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

typedef struct riverbc_struct
{
    double          head;
    double          flux;
} riverbc_struct;

typedef struct riveric_struct
{
    double          stage;
    double          gw;
} riveric_struct;

typedef struct river_struct
{
    attrib_struct   attrib;
    river_topo_struct     topo;
    shp_struct      shp;
    matl_struct     matl;
    river_wstate_struct   ws;
    river_wstate_struct   ws0;
    river_wflux_struct    wf;
    riveric_struct  ic;
    riverbc_struct  bc;
#ifdef _DAILY_
    river_daily_struct    daily;
#endif
    int             leftele;    /* Left neighboring element */
    int             rightele;   /* Right neighboring element */
    int             fromnode;   /* Upstream Node no. */
    int             tonode;     /* Dnstream Node no. */
    int             up;
    int             down;       /* down stream segment */
#ifdef _CYCLES_
    solute_struct   NO3sol;
    solute_struct   NH4sol;
#endif
#ifdef _BGC_
    river_stor_struct     stor;     /* meteorological data array */
    double          sminn;
    double          nleached_snk;
    double          sminn_leached;
#endif
} river_struct;

#endif
