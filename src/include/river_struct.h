#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/* River attribute */
typedef struct river_attrib_struct
{
    int             riverbc_type;           /* river boundary condition type */
} river_attrib_struct;

/* River topography parameters */
typedef struct river_topo_struct
{
    double          area;                   /* area (m2) */
    double          x;                      /* x of centroid (m) */
    double          y;                      /* y of centroid (m) */
    double          zmax;                   /* river bank elevation (m) */
    double          zbed;                   /* river bed elevation (m) */
    double          node_zmax;              /* elevation of the downstream node
                                             * (m) */
    double          dist_left;              /* distance to left neighbor (m) */
    double          dist_right;             /* distance to right neighbor (m) */
} river_topo_struct;

/* River water states */
typedef struct river_wstate_struct
{
    double          stage;                  /* river stage (m) */
} river_wstate_struct;

/* River water fluxes */
typedef struct river_wflux_struct
{
    double          rivflow[NUM_RIVFLX];    /* river fluxes (m3 s-1) */
} river_wflux_struct;

/* River shape parameters */
typedef struct shp_struct
{
    double          depth;                  /* river channel depth (m) */
    int             intrpl_ord;             /* interpolation order (shape of
                                             * channel) */
    double          coeff;                  /* width coefficient */
    double          length;                 /* length of channel (m) */
    double          width;                  /* width of channel (m) */
} shp_struct;

typedef struct matl_struct
{
    double          rough;                  /* river channel roughness (s m-1/3)
                                             */
    double          cwr;                    /* discharge coefficient (-) */
    double          ksath;                  /* bank hydraulic conductivity
                                             * (m s-1) */
#if defined(_CYCLES_OBSOLETE_)
    double          bd;
#endif
} matl_struct;

/* River boundary conditions */
typedef union river_bc_struct
{
    union
    {
        double          head;               /* value of Dirichlet-type boundary
                                             * condition (m) */
        double          flux;               /* value of Neumann-type boundary
                                             * condition (m3 s-1) */
    };
} river_bc_struct;

/* River initial conditions */
typedef struct river_ic_struct
{
    double          stage;
} river_ic_struct;

#if defined(_BGC_) && !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
/* River nitrogen state variables */
typedef struct river_nstate_struct
{
    double          streamn;                /* stream N pool (kgN m-2) */
} river_nstate_struct;

/* Daily river nitrogen flux variables */
typedef struct river_nflux_struct
{
    double          sminn_leached;          /* leaching flux (kgN m-2 day-1) */
} river_nflux_struct;

/* River CN initial conditions */
typedef struct river_bgcic_struct
{
    double          streamn;
} river_bgcic_struct;
#endif

#if defined(_CYCLES_)
typedef struct river_nstate_struct
{
    double          no3;
    double          nh4;
} river_nstate_struct;
#endif

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
typedef struct river_solute_struct
{
    double          conc;                   /* solute concentration */
    double          flux[NUM_RIVFLX];       /* solute flux (mass or mol s-1) */
} river_solute_struct;
#endif

/* River structure */
typedef struct river_struct
{
    int             ind;                    /* river index */
    int             leftele;                /* left neighbor*/
    int             rightele;               /* right neighbor */
    int             fromnode;               /* upstream node */
    int             tonode;                 /* downstream node */
    int             down;                   /* down stream channel segment */
    river_attrib_struct attrib;
    river_topo_struct topo;
    shp_struct      shp;
    matl_struct     matl;
    river_wstate_struct ws;
    river_wstate_struct ws0;
    river_wflux_struct wf;
    river_ic_struct ic;
    river_bc_struct bc;
#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
    river_solute_struct solute[NSOLUTE];
#endif
#if defined(_CYCLES_)
    river_nstate_struct ns;
#endif
#if defined(_BGC_) && !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
    river_nstate_struct ns;
    river_nflux_struct nf;
    river_solute_struct nsol;
    river_bgcic_struct restart_input;
    river_bgcic_struct restart_output;
#endif
#if defined(_RT_)
    chmstate_struct     chms;
#endif
} river_struct;

#endif
