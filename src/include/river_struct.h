#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/*****************************************************************************
 * River attribute
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * riverbc_type             int         river boundary condition type
 ****************************************************************************/
typedef struct river_attrib_struct
{
    int             riverbc_type;
} river_attrib_struct;

/*****************************************************************************
 * River topography parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * area                     double      area of element [m2]
 * x                        double      x of centroid [m]
 * y                        double      y of centroid [m]
 * zmin                     double      bedrock elevation [m]
 * zmax                     double      river bank elevation [m]
 * zbed                     double      river bed elevation [m]
 * node_zmax                double      elevation of the downstream node [m]
 * dist_left                double      distance to left neighbor [m]
 * dist_right                double      distance to right neighbor [m]
 ****************************************************************************/
typedef struct river_topo_struct
{
    double          area;
    double          x;
    double          y;
    double          zmin;
    double          zmax;
    double          zbed;
    double          node_zmax;
    double          dist_left;
    double          dist_right;
} river_topo_struct;

/*****************************************************************************
 * River water states
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * stage                    double      river stage [m]
 * gw                       double      groundwater level [m]
 ****************************************************************************/
typedef struct river_wstate_struct
{
    double          stage;
    double          gw;
} river_wstate_struct;

/*****************************************************************************
 * River water fluxes
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * rivflow                  double[]    river fluxes [m3 s-1]
 ****************************************************************************/
typedef struct river_wflux_struct
{
    double          rivflow[NUM_RIVFLX];
} river_wflux_struct;

/*****************************************************************************
 * River shape parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * depth                    double      river channel depth [m]
 * intrpl_ord               int         interpolation order (shape of channel)
 * shpcoeff                 double      width coefficient
 * length                   double      length of channel [m]
 * width                    double      width of channel [m]
 ****************************************************************************/
typedef struct shp_struct
{
    double          depth;
    int             intrpl_ord;
    double          coeff;
    double          length;
    double          width;
} shp_struct;

/*****************************************************************************
 * River channel matierial parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * rough                    double      river channel roughness [s m-1/3]
 * cwr                      double      discharge coefficient [-]
 * ksath                    double      bank hydraulic conductivity [m s-1]
 * ksatv                    double      bed hydraulic conductivity [m s-1]
 * bedthick                 double      bed thickness [m]
 * porosity                 double      bed porosity [m3 m-3]
 ****************************************************************************/
typedef struct matl_struct
{
    double          rough;
    double          cwr;
    double          ksath;
    double          ksatv;
    double          bedthick;
    double          porosity;
} matl_struct;

/*****************************************************************************
 * River boundary conditions
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * head                     double      value of Dirichlet-type boudnary
 *                                        condition [m]
 * flux                     double      value of Neumann-type boundary
 *                                        condition [m3 s-1]
 ****************************************************************************/
typedef struct river_bc_struct
{
    double          head;
    double          flux;
} river_bc_struct;

/*****************************************************************************
 * River initial conditions
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * stage                    double      river stage [m]
 * gw                       double      groundwater level [m]
 ****************************************************************************/
typedef struct river_ic_struct
{
    double          stage;
    double          gw;
} river_ic_struct;

#ifdef _BGC_
/*****************************************************************************
 * River nitrogen state variables
 * ---------------------------------------------------------------------------
 * streamn                  double      stream N [kgN m-2]
 * sminn                    double      soil mineral N [kgN m-2]
 ****************************************************************************/
typedef struct river_nstate_struct
{
    double          streamn;
    double          sminn;
} river_nstate_struct;

/*****************************************************************************
 * Daily river nitrogen flux variables
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * sminn_leached            double      leaching flux (kgN/m2/d)
 ****************************************************************************/
typedef struct river_nflux_struct
{
    double          sminn_leached;
} river_nflux_struct;

typedef struct river_solute_struct
{
    double          conc_stream;
    double          conc_bed;
    double          flux[NUM_RIVFLX];
} river_solute_struct;

typedef struct river_bgcic_struct
{
    double          streamn;
    double          sminn;
} river_bgcic_struct;

#endif

#ifdef _CYCLES_
typedef struct river_solute_struct
{
    double          soluteMass;
    double          soluteMassAdsorbed;
    double          soluteConc;
    double          soluteFluxLat[4];
} river_solute_struct;
#endif

/*****************************************************************************
 * River structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * leftele                  int         Left bank id
 * rightele                 int         right bank id
 * fromnode                 int         upstream node id
 * tonode                   int         downstream node id
 * up                       int         upstream channel id
 * down                     int         downstream channel id
 ****************************************************************************/
typedef struct river_struct
{
    int             leftele;    /* Left neighboring element */
    int             rightele;   /* Right neighboring element */
    int             fromnode;   /* Upstream Node no. */
    int             tonode;     /* Dnstream Node no. */
    int             up;
    int             down;       /* down stream segment */
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
    river_solute_struct NO3sol;
    river_solute_struct NH4sol;
    double          NO3Leaching[4];
    double          NH4Leaching[4];
#endif
#ifdef _BGC_
    river_nstate_struct ns;
    river_nflux_struct nf;
    river_solute_struct nsol;
    river_bgcic_struct restart_input;
    river_bgcic_struct restart_output;
#endif
} river_struct;
#endif
