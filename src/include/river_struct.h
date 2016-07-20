#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/* Attribute */
typedef struct river_attrib_struct
{
    int             riverbc_type;
} river_attrib_struct;

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
    double          rivflow[NUM_RIVFLX];
} river_wflux_struct;

#ifdef _BGC_
typedef struct river_stor_struct
{
    double         *stage;
    double         *gw;
    double         *rivflow[NUM_RIVFLX];
    int            *flag;
} river_stor_struct;
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

typedef struct river_bc_struct
{
    double          head;
    double          flux;
} river_bc_struct;

typedef struct river_ic_struct
{
    double          stage;
    double          gw;
} river_ic_struct;

#ifdef _DAILY_
typedef struct river_daily_struct
{
    int             counter;

    double          avg_stage;
    double          avg_gw;

    double          avg_rivflow[NUM_RIVFLX];
} river_daily_struct;
#endif

#ifdef _BGC_
typedef struct river_nstate_struct
{
    double          sminn;
} river_nstate_struct;

typedef struct river_nflux_struct
{
    double          sminn_leached;
} river_nflux_struct;
#endif

typedef struct river_struct
{
    river_attrib_struct attrib;
    river_topo_struct topo;
    shp_struct      shp;
    matl_struct     matl;
    river_wstate_struct ws;
    river_wstate_struct ws0;
    river_wflux_struct wf;
    river_ic_struct ic;
    river_bc_struct bc;
#ifdef _DAILY_
    river_daily_struct daily;
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
    river_stor_struct stor;     /* meteorological data array */
    river_nstate_struct ns;
    river_nflux_struct nf;
#endif
} river_struct;
#endif
