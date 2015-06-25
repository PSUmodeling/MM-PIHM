#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

/* Data model for a triangular element */
typedef struct element_type
{
    int             index;      /* Element No. */
    int             node[3];    /* Counterclock-wise */
    int             nabr[3];    /* neighbor i shares edge i
                                 * (0: on boundary) */

    realtype        edge[3];    /* edge i is from node i to node i+1 */
    realtype        area;       /* area of element */

    realtype        x;          /* x of centroid */
    realtype        y;          /* y of centroid */
    realtype        zmin;       /* z_min of centroid */
    realtype        zmax;       /* z_max of centroid */

    realtype        KsatH;      /* horizontal geologic saturated
                                 * hydraulic conductivity */
    realtype        KsatV;      /* vertical geologic saturated
                                 * hydraulic conductivity */
    realtype        infKsatV;   /* vertical surface saturated hydraulic
                                 * conductivity */
    realtype        Porosity;
    realtype        infD;       /* depth from ground surface accross
                                 * which head is calculated during
                                 * infiltration */
    realtype        Alpha;      /* Alpha from van Genuchten eqn */
    realtype        Beta;
    realtype        ThetaS;
    realtype        ThetaR;
    realtype        ThetaRef;   /* YS: Soil field capacity */
    realtype        ThetaW;     /* YS: Soil wilting point */
    realtype        RzD;        /* Root zone depth */
    realtype        macD;       /* macropore Depth */
    realtype        macKsatH;   /* macropore horizontal saturated
                                 * hydraulic conductivity */
    realtype        macKsatV;   /* macropore vertical saturated
                                 * hydraulic conductivity */
    realtype        vAreaF;     /* macropore area fraction on a
                                 * vertical cross-section */
    realtype        hAreaF;     /* macropore area fraction on a
                                 * horizontal cross-section */
    int             Macropore;  /* 1: macropore; 0: regular soil */

    realtype        LAImax;     /* maxm. LAI accross all seasons for a
                                 * vegetation type */
    realtype        VegFrac;    /* areal vegetation fraction in a
                                 * triangular element */
    realtype        Rs_ref;     /* reference incoming solar flux for
                                 * photosynthetically active canopy */
    realtype        Rmin;       /* minimum canopy resistance */
    realtype        Rough;      /* surface roughness of an element */

    realtype        windH;      /* wind measurement height */
#ifdef _RT_
    realtype        temp;       /* temperature   */
#endif
    int             soil;       /* soil type */
    int             geol;       /* geology type */
    int             LC;         /* Land Cover type  */
    int             IC;         /* initial condition type */
    int             BC[3];      /* Boundary condition type.
                                 * 0: Natural BC (no flow);
                                 * 1: Dirichlet BC;
                                 * 2:Neumann BC */
    int             meteo;      /* precipitation (forcing) type */
    int             LAI;        /* LAI forcing type (0: use climatological
                                 * values; else: use user provided time
                                 * series */
    int             source;     /* source (well) type */
    int             meltF;      /* meltFactor */

    /*
     * for calculation of dh/ds
     */
    realtype        surfH[3];   /* Total head in neighboring cells */
    realtype        surfX[3];   /* Center X location of neighboring
                                 * cells */
    realtype        surfY[3];   /* Center Y location of neighboring
                                 * cells */
    realtype        dhBYdx;     /* Head gradient in x dirn. */
    realtype        dhBYdy;     /* Head gradient in y dirn. */
} element;

/*
 * Data model for a node
 */
typedef struct nodes_type
{
    int             index;      /* Node no. */
    realtype        x;          /* x coordinate */
    realtype        y;          /* y coordinate */
    realtype        zmin;       /* z bed rock elevation */
    realtype        zmax;       /* z surface elevation  */
} nodes;

/*
 * Initial state variable conditions on each element
 */
typedef struct element_IC_type
{
    int             index;
    realtype        interception;   /* Interception storage (Note all these
                                     * variables have dimension of L */
    realtype        snow;       /* Snow depth */
    realtype        surf;       /* Overland flow depth */
    realtype        unsat;      /* unsaturated zone depth */
    realtype        sat;        /* saturated zone depth */
} element_IC;

typedef struct soils_type
{
    int             index;      /* index */
    realtype        KsatV;      /* vertical saturated soil
                                 * conductivity */
    realtype        ThetaS;     /* soil porosity */
    realtype        ThetaR;     /* soil moisture residual */
    realtype        ThetaW;     /* YS: wilting point */
    realtype        ThetaRef;   /* YS: field capacity */
    realtype        qtz;        /* YS: quartz content */
    realtype        Alpha;      /* soil curve parameter 1 */
    realtype        Beta;       /* soil curve parameter 2 */

    realtype        hAreaF;     /* macroporous area fraction on
                                 * horizontal section */
    realtype        macKsatV;   /* macroporous saturated vertical
                                 * conductivity */
    realtype        infD;       /* depth from ground surface accross which
                                 * head is calculated during infiltration */
} soils;

typedef struct geol_type
{
    int             index;      /* index */
    realtype        KsatH;      /* horizontal saturated geology
                                 * conductivity */
    realtype        KsatV;      /* vertical saturated geology
                                 * conductivity */
    realtype        ThetaS;     /* geology porosity */
    realtype        ThetaR;     /* residual porosity */
    realtype        Alpha;      /* van Genuchten parameter */
    realtype        Beta;       /* van Genuchten parameter */

    realtype        vAreaF;     /* macroporous area fraction on vertical
                                 * section */
    realtype        macKsatH;   /* macroporous saturated
                                 * horizontal conductivity */
    realtype        macD;
} geol;

typedef struct lc_type
{
    int             index;      /* index */

    realtype        LAImax;     /* max LAI */
    realtype        LAImin;     /* YS: min LAI */
    realtype        VegFrac;    /* Canopy Fracn */
    realtype        Albedo;     /* Albedo */
    realtype        Albedo_min; /* YS: Minimum albedo */
    realtype        Albedo_max; /* YS: Maximum albedo */
    realtype        Emiss_min;  /* YS: Minimum emissivity */
    realtype        Emiss_max;  /* YS: Maximum emissivity */
    realtype        z0_min;     /* YS: Minimum roughness length */
    realtype        z0_max;     /* YS: Maximum roughness length */
    realtype        h_s;        /* YS: Vapor pressure deficit stress
                                 * parameter */
    realtype        snup;       /* YS */
    realtype        Rs_ref;     /* Visible solar flux used in radiation
                                 * stress */
    realtype        Rmin;       /* Minimum stomatal resistance */
    realtype        Rough;      /* Surface roughness factor  */
    realtype        RzD;        /* RootZone Depth */
} LC;

typedef struct river_segment_type
{
    int             index;
    realtype        x;          /* x of river segment */
    realtype        y;
    realtype        zmin;       /* bed elevation  */
    realtype        zmax;       /* bank elevation */
    realtype        depth;      /* max depth */
    realtype        Length;     /* Riv segment Length */
    realtype        Rough;      /* Manning's roughness coeff */
    realtype        KsatH;      /* Side conductivity */
    realtype        KsatV;      /* Bed conductivity */
    realtype        bedThick;
    realtype        coeff;      /* Coefficient c in
                                 * D=c*pow(B/2, interpOrd),where D is depth */
    int             FromNode;   /* Upstream Node no. */
    int             ToNode;     /* Dnstream Node no. */
    int             down;       /* down stream segment */
    int             LeftEle;    /* Left neighboring element */
    int             RightEle;   /* Right neighboring element */
    int             shape;      /* shape type    */
    int             material;   /* material type */
    int             IC;         /* IC type */
    int             BC;         /* BC type */
    int             reservoir;

} river_segment;

typedef struct river_shape_type
{
    int             index;
    realtype        depth;      /* depth */
    int             interpOrd;  /* Interpolation order for river shape:
                                 * 1: rectangle,
                                 * 2: triangle,
                                 * 3: quadratic,
                                 * 4: cubic */
    realtype        coeff;      /* Coefficient c in
                                 * D = c * pow(B / 2, interpOrd) */
} river_shape;

typedef struct river_material_type
{
    int             index;
    realtype        Rough;
    realtype        Cwr;        /* Weir Discharge Coefficient */
    realtype        KsatH;      /* Conductivity of river banks */
    realtype        KsatV;      /* Conductivity of river bed */
    realtype        bedThick;   /* thickeness of conductive river bed */
} river_material;

typedef struct river_IC_type
{
    int             index;
    realtype        value;      /* initial flow depth */
} river_IC;


/* Global calibration sturcture
 * For explanation of each calibration variable, look for corresponding
 * variables above */
typedef struct global_calib
{
    realtype        KsatH;
    realtype        KsatV;
    realtype        infKsatV;
    realtype        macKsatH;
    realtype        macKsatV;
    realtype        infD;
    realtype        RzD;
    realtype        macD;
    realtype        Porosity;
    realtype        Alpha;
    realtype        Beta;
    realtype        vAreaF;
    realtype        hAreaF;
    realtype        Temp;
    realtype        Prep;
    realtype        VegFrac;
    realtype        Albedo;
    realtype        Rough;

    realtype        rivRough;
    realtype        rivKsatH;
    realtype        rivKsatV;
    realtype        rivbedThick;
    realtype        rivDepth;
    realtype        rivShapeCoeff;

    realtype        ThetaRef;   /* YS */
    realtype        ThetaW;     /* YS */
    realtype        Rmin;       /* YS */
#ifdef _FLUX_PIHM_
    realtype        TF;         /* YS */
    realtype        IS;         /* YS */
    realtype        Czil;       /* YS */
    realtype        fx_soil;    /* YS */
    realtype        fx_canopy;  /* YS */
    realtype        Rs_ref;     /* YS */
    realtype        h_s;        /* YS */
#endif
#ifdef _RT_
  realtype   PCO2;
  realtype   Keq;
  realtype   SSA;
  realtype   Site_den;
  realtype   Prep_conc;
#endif
} globalCal;

typedef struct process_control
{
    realtype        Et0;
    realtype        Et1;
    realtype        Et2;
} processCal;

/* Model_data definition */
typedef struct model_data_structure
{
    int             UnsatMode;  /* Unsat Mode */
    int             SurfMode;   /* Surface Overland Flow Mode */
    int             RivMode;    /* River Routing Mode */

    int             NumEle;     /* Number of Elements */
    int             NumNode;    /* Number of Nodes */
    int             NumRiv;     /* Number of Rivere Segments */

    int             NumLAI;
    int             NumSource;  /* Number of Source time series types */
    int             NumMeltF;   /* Number of Melt Factor Time series */

    int             NumSoil;    /* Number of Soils */
    int             NumGeol;    /* Number of Geologies */
    int             NumRes;     /* Number of Reservoir */
    int             NumLC;      /* Number of Land Cover Index Data */


    int             Num1BC;     /* Number of Dirichlet BC */
    int             Num2BC;     /* Number of Numann BC */
    int             NumBC;
    int             NumEleIC;   /* Number of Element Initial Condtion */

    int             NumRivShape;    /* Number of River Shape */
    int             NumRivMaterial; /* Number of River Bank/Bed Material */
    int             NumRivIC;   /* Number of River Initial Condition */
    int             NumRivBC;   /* Number of River Boundary Condition */

    realtype        Rmax;       /* YS */
    int             bare;       /* YS */
    realtype        fx_canopy;  /* YS */
    realtype        Tref;       /* YS */

    element        *Ele;        /* Store Element Information */
    nodes          *Node;       /* Store Node Information */
    element_IC     *Ele_IC;     /* Store Element Initial Condtion */
    soils          *Soil;       /* Store Soil Information */
    geol           *Geol;       /* Store Geology Information */
    LC             *LandC;      /* Store Land Cover Information */

    river_segment  *Riv;        /* Store River Segment Information */
    river_shape    *Riv_Shape;  /* Store River Shape Information */
    river_material *Riv_Mat;    /* Store River Bank Material Information */
    river_IC       *Riv_IC;     /* Store River Initial Condition */

    realtype       *ISFactor;   /* ISFactor is used to calculate ISMax from
                                 * LAI */
    realtype       *windH;      /* Height at which wind speed is measured */

    TSD            *TSD_EleBC;  /* Element Boundary Condition Time Series
                                 * Data */
    TSD            *TSD_meteo;    /* YS */
    TSD            *TSD_lai;    /* YS */
    TSD            *TSD_rl;
    TSD            *TSD_mf;     /* YS  */
    TSD            *TSD_ss;     /* YS */
    int             NumTS;      /* YS */
    TSD            *TSD_Riv;    /* River Related Time Series Data  */

    realtype      **FluxSurf;   /* Overland Flux */
    realtype      **FluxSub;    /* Subsurface Flux */
    realtype      **FluxRiv;    /* River Segement Flux */
#ifdef _RT_
    realtype      **AreaSub;    /* Area of contact between saturated cells */
#endif
 
    realtype       *ElePrep;    /* Precep. on each element */
    realtype       *EleNetPrep; /* Net precep. on each elment */
    realtype       *EleViR;     /* Variable infiltration rate */
    realtype       *Recharge;   /* Recharge rate to GW */
    realtype       *EleSnow;    /* YS: Snow water equivalent on each element*/
    realtype       *EleSnowGrnd;    /* Snow depth on ground element */
    realtype       *EleSnowCanopy;  /* Snow depth on canopy element */
    realtype       *EleIS;      /* Interception storage */
    realtype       *EleISmax;   /* Maximum interception storage
                                 * (liquid precep) */
    realtype       *EleISsnowmax;   /* Maximum snow interception storage */
    realtype       *EleTF;      /* Through Fall */
    realtype      **EleET;      /* Evapotranspiration (from canopy, ground,
                                 * transpiration) */

    realtype       *Albedo;     /* albedo of a triangular element */

    realtype       *EleSurf;    /* YS: Stores surface water level of
                                 * last time step */
    realtype       *RivStg;     /* YS: Stores river stage of last time step */
    realtype       *EleGW;      /* YS: Stores groundwater level of last time
                                 * step */
    realtype       *EleUnsat;   /* YS: Stores unsaturated storage of
                                 * last time step */
    int            *EleMacAct;
#ifdef _FLUX_PIHM_
    realtype       *SfcSat;     /* YS: Surface saturation */
    realtype       *EleETsat;   /* YS: Fraction of Transpiration that
                                 * extracts from the saturated zone */
    realtype       *EleFCR;     /* YS: reduction of infiltration caused
                                 * by frozen ground */
    realtype	   *avg_inf;	/* YS: average infiltration over an ET time
				 * interval */ 
    realtype	   *avg_rech;	/* YS: average recharge over an ET time interval */
    realtype	  **avg_subflux;/* YS: average subsurfae flux over an ET time
				     * interval */
#endif

#ifdef _BGC_
    realtype       *daily_tmax;
    realtype       *daily_tmin;
    realtype       *daily_tnight;
    realtype       *daily_tavg;
    realtype       *daily_tday;
    realtype       *daily_tsoil;
    realtype       *daily_swc;
    realtype       *daily_soilw;
    realtype       *daily_q2d;
    realtype       *daily_swavgfd;
    realtype       *daily_par;
    realtype       *daily_pa;
    realtype      **daily_subflux;
    realtype       *daily_albedo;
    realtype       *daily_ch;
#endif

    realtype       *DummyY;
    processCal      pcCal;

    realtype        dt;         /* YS: Time step */
} *Model_Data;

typedef struct control_data_structure
{
    int             Verbose;
    int             Debug;
    int             Ascii;      /* YS: Add Ascii output model
                                 * (default is binary */
    int             Spinup;     /* YS: Runs model as spinup. Model output at
                                 * the last step will be saved in .init */
    int             Solver;     /* Solver type */
    int             NumSteps;   /* Number of external time steps (when
                                 * results can be printed) for the
                                 * whole simulation */
    int             NumPrint;   /* YS: Number of variables for output */

    /* Time interval to output average values of variables
     * Variables will not be printed if time intervals are set to 0 */
    int             PrintGW;
    int             PrintSurf;
    int             PrintSnow;
    int             PrintRivStg;
    int             PrintRech;
    int             PrintIS;
    int             PrintUnsat;
    int             PrintET[3];
    int             PrintRivFlx[11];
    int             PrintSubFlx;

    Print_Ctrl      PCtrl[100]; /* YS */

    int             init_type;  /* initialization mode */

    realtype        abstol;     /* absolute tolerance */
    realtype        reltol;     /* relative tolerance */
    realtype        InitStep;   /* initial step size */
    realtype        MaxStep;    /* Maximum step size */
    realtype        ETStep;     /* Step for et from interception */

    int             GSType;
    int             MaxK;       /* Maximum Krylov order */
    realtype        delt;

    realtype        StartTime;  /* Start time of simulation */
    realtype        EndTime;    /* End time of simulation */

    int             outtype;
    realtype        a;          /* External time stepping controls */
    realtype        b;

    realtype       *Tout;

    globalCal       Cal;        /* Convert this to pointer for
                                 * localized calibration */
} * Control_Data;
#endif
