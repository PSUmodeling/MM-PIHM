#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/* River attribute */
typedef struct river_attrib_struct
{
    int             riverbc_type;    /* river boundary condition type */
} river_attrib_struct;

/* River topography parameters */
typedef struct river_topo_struct
{
    double          area;          /* area of element (m2) */
    double          x;             /* x of centroid (m) */
    double          y;             /* y of centroid (m) */
    double          zmin;          /* bedrock elevation (m) */
    double          zmax;          /* river bank elevation (m) */
    double          zbed;          /* river bed elevation (m) */
    double          node_zmax;     /* elevation of the downstream node (m) */
    double          dist_left;     /* distance to left neighbor (m) */
    double          dist_right;    /* distance to right neighbor (m) */
} river_topo_struct;

/* River water states */
typedef struct river_wstate_struct
{
    double          stage;    /* river stage (m) */
    double          gw;       /* groundwater level (m) */
} river_wstate_struct;

/* River water fluxes */
typedef struct river_wflux_struct
{
    double          rivflow[NUM_RIVFLX];    /* river fluxes (m3 s-1) */
} river_wflux_struct;

/* River shape parameters */
typedef struct shp_struct
{
    double          depth;         /* river channel depth (m) */
    int             intrpl_ord;    /* interpolation order (shape of channel) */
    double          coeff;         /* width coefficient */
    double          length;        /* length of channel (m) */
    double          width;         /* width of channel (m) */
} shp_struct;

typedef struct matl_struct
{
    double          rough;       /* river channel roughness (s m-1/3) */
    double          cwr;         /* discharge coefficient (-) */
    double          ksath;       /* bank hydraulic conductivity (m s-1) */
    double          ksatv;       /* bed hydraulic conductivity (m s-1) */
    double          bedthick;    /* bed thickness (m) */
    double          porosity;    /* bed porosity (m3 m-3) */
    double          smcmin;      /* bed residual soil moisture content (m3 m-3)
                                  */
} matl_struct;

/* River boundary conditions */
typedef struct river_bc_struct
{
    double          head;   /* value of Dirichlet-type boundary condition (m) */
    double          flux;   /* value of Neumann-type boundary condition (m3 s-1)
                             */
} river_bc_struct;

/* River initial conditions */
typedef struct river_ic_struct
{
    double          stage;
    double          gw;
} river_ic_struct;

#if defined(_BGC_) && !defined(_LUMPED_)
/* River nitrogen state variables */
typedef struct river_nstate_struct
{
    double          streamn;    /* stream N pool (kgN m-2) */
    double          sminn;      /* river bed soil mineral N (kgN m-2) */
} river_nstate_struct;

/* Daily river nitrogen flux variables */
typedef struct river_nflux_struct
{
    double          sminn_leached;    /* leaching flux (kgN m-2 day-1) */
} river_nflux_struct;

/* River solute transport structure */
typedef struct river_solute_struct
{
    double          conc_stream;         /* stream pool concentration
                                          * (kg kgH2O-1) */
    double          conc_bed;            /* bed pool concentration (kg kgH2O-1)
                                          */
    double          flux[NUM_RIVFLX];    /* solute fluxes (kg s-1) */
} river_solute_struct;

/* River CN initial conditions */
typedef struct river_bgcic_struct
{
    double          streamn;
    double          sminn;
} river_bgcic_struct;
#endif

#ifdef _CYCLES_
/* River daily average */
typedef struct river_daily_struct
{
    int             counter;                /* counter used for averaging */
    double          avg_stage;              /* daily average river stage (m) */
    double          avg_gw;                 /* daily average groundwater level
                                             * (m) */
    double          avg_rivflow[NUM_RIVFLX];/* daily average river flux (m3 s-1)
                                             */
} river_daily_struct;

typedef struct river_solute_struct
{
    double          soluteMass;
    double          soluteMassAdsorbed;
    double          soluteConc;
    double          soluteFluxLat[4];
} river_solute_struct;

typedef struct river_cyclesic_struct
{
    double          NO3_Mass;
    double          NH4_Mass;
} river_cyclesic_struct;
#endif

/* River structure */
typedef struct river_struct
{
    int             ind;         /* river index */
    int             leftele;     /* left neighboring element */
    int             rightele;    /* right neighboring element */
    int             fromnode;    /* upstream node */
    int             tonode;      /* downstream node */
    int             down;        /* down stream channel segment */
    river_attrib_struct attrib;
    river_topo_struct topo;
    shp_struct      shp;
    matl_struct     matl;
    river_wstate_struct ws;
    river_wstate_struct ws0;
    river_wflux_struct wf;
    river_ic_struct ic;
    river_bc_struct bc;
#ifdef _CYCLES_
    river_daily_struct daily;
    river_solute_struct NO3sol;
    river_solute_struct NH4sol;
    river_cyclesic_struct cycles_restart;
    double          NO3Leaching[4];
    double          NH4Leaching[4];
#endif
#if defined(_BGC_) && !defined(_LUMPED_)
    river_nstate_struct ns;
    river_nflux_struct nf;
    river_solute_struct nsol;
    river_bgcic_struct restart_input;
    river_bgcic_struct restart_output;
#endif
} river_struct;
#endif
