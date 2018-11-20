/******************************************************************************
 * File        : rt.h
 * Function    : Declaration and Definition of reaction variables & data struc\
 *               ture. Data could be mainly be classified into three classes: \
 *               chemical knowledge, concentration distribution and grid struc\
 *               ture.
 * Developer of PIHM RT 1.0: Chen Bao (baochen.d.s@gmail.com)
 * Date        : June 2013
 *****************************************************************************/

/* parameters of Debye_Huckel Equation */
typedef struct Debye_Huckel_structure
{
    double          adh;
    double          bdh;
    double          bdt;
} Debye_Huckel;

typedef struct vol_conc_type
{
    int             index;      /* Volume No. Note this number may be different than the element No. in PIHM */
    int             illness;    /* black list */
    int             NumStc;     /* Total Number of primary (total) species */
    int             NumSsc;     /* Total Number of secondary species */
    int             type;       /* type of volume: unsaturated, groundwater, river, or river bed */
    double         *t_conc;     /* Concentration of species x, default unit in Mol/kg water. */
    double         *p_conc;     /* Primary concentrations, default unit in Mol/kg water */
    /* It should be noted that this array does not only store aqueous concentration
     * but also include other types of primary species that are not mobile
     * These various types of primary species are strictly stored in certain order so that
     * enables easier implementation of OS3D scheme */
    /* other units:
     * for surface complexation
     * for cation exchange
     * for minerals
     */
    double         *s_conc;     /* Secondary concentrations, default unit in Mol/kg water */
    double         *p_actv;     /* Activity of primary species */
    double         *p_para;     /* parameters of primary species
                                 * for surface complexation
                                 * for cation exchage
                                 * for minerals
                                 */
    int            *p_type;     /* type of primary species,
                                 * 1: primary aqueous
                                 * 2: primary adsorption
                                 * 3: primary cation exchange
                                 * 4: primary mineral
                                 */
    double          sat;        /* Saturation of Volume, dimensionless, from 0 to 1 */
    double          height_o;   /* height of volume, at previous time step */
    double          height_t;   /* height of volume, at current time step */
    double          height_int; /* temporary variable for intrapolation of gw height */
    double          height_sp;  /* slope of the height change during this period */
    double          height_v;   /* height of the total block, static */
    double          maxwater;   /* maximum water height during 2016, by Wei 01.21 */
    double          area;       /* area of the triangular element */
    double          vol_o;      /* volume of water of last time step */
    double          vol;        /* volume of water of current time step */
    double          porosity;   /* porosity of the volume */
    double          rt_step;    /* rt_step of cell */
    double         *log10_pconc;     /* for output only */
    double         *log10_sconc;     /* for output only */
    double         *btcv_pconc; /* for btcv output only */
} vol_conc;

typedef struct face_type
{
    int             nodeup;     /* Volume No. of first node */
    int             nodelo;     /* Volume No. of second node */
    int             nodeuu;     /* Volume No. of node before first node */
    int             nodell;     /* Volume No. of node after second node */
    int             node_trib;  // # node of tributary; > 0 indicates a tributary inflow, 01.14 by Wei Zhi
    int             BC;         /* This face is a boundary face, do not do dispersion and diffusion */
    double          distance;   /* distance from centroid of first node to second node, in meter */
    double          velocity;   /* linear velocity of flux at this face, in m/d */
    double          flux;       /* flux at surface, in cubic m per day  */
    double          flux_trib;  // tributary fluc, 01.14 by Wei Zhi
    double          s_area;     /* contact surface area, in square m */
} face;

typedef struct species_type
{
    double          DiffCoe;    /* diffusion coefficient, measured in cm2/s */
    double          DispCoe;    /* dispersion coefficient, measured in cm2/s */
    double          MolarMass;  /* measured in g/mol */
    double          MolarVolume;    /* measured in cm3/mol */
    double          Charge;     /* Array of species charge */
    double          SizeF;      /* Array of species size factor for DH equation */
    char           *ChemName;   /* usually the molecular formula of the very chemical, or any convenient name as long as the database has it */
    int             itype;      /* type of primary species,
                                 * 1: primary aqueous
                                 * 2: primary adsorption
                                 * 3: primary cation exchange
                                 * 4: primary mineral
                                 */
    int             mtype;      /* type of the mass action species
                                 * 1: mobile mass action species
                                 * 0: immobile mass action species
                                 * 2: mixed mobility mass action species
                                 */
} species;

typedef struct Kinetic_Reaction_structure
{
    char            species[MAXSTRING];     // target mineral
    int             position;               // position of target mineral in the
                                            // array of primary species
    char            Label[MAXSTRING];       // label of kinetic reaction, could
                                            // be multiple for one mineral
                                            // (aqueous species)
    int             type;                   // type of the kinetic reaction
                                            // 1: tst  2: precipitation-only
                                            // 3: dissolution-only 4: monod
    double          rate;                   // rate of kinetic reaction, often
                                            // at 25 C
    double          actv;                   // activation energy, used to
                                            // calculate kinetic rate of
                                            // reaction under different
                                            // temperatures
    double          Keq;                    // equilibrium constant
    int             num_dep;                // number of dependency
    char            dep_species[MAXDEP][MAXSTRING]; // species that this kinetic
                                            // rate depends on
    int             dep_position[MAXDEP];   // store position of species that
                                            // this kinetic reaction depends on
    double          dep_power[MAXDEP];      // power of dependency
    char            biomass_species[MAXSTRING];
    int             biomass_position;
    int             num_monod;
    char            monod_species[MAXDEP][MAXSTRING]; // species that this
                                            // kinetic rate depends on for monod
                                            // type
    int             monod_position[MAXDEP];
    double          monod_para[MAXDEP];     // parameters for monod dependency
    int             num_inhib;
    char            inhib_species[MAXDEP][MAXSTRING]; // species that this
                                            // kinetic rate is inhibited by, for
                                            // monod type
    int             inhib_position[MAXDEP];
    double          inhib_para[MAXDEP];     // parameters that controls this
                                            // inhibition
} Kinetic_Reaction;


typedef struct Pump_Data_structure
{
    int             Pump_Location;  /* Index of pump grid block */
    int             Position_Species;   /* Index of chemical species */
    char           *Name_Species;   /* The name of chemical species */
    double          Injection_rate; /* Rate of injection, or extraction, in moles/year */
    double          Injection_conc; /* Concentration of injection in moles/L */
    double          flow_rate;  /* Calculated flow rate from the above two parameters */
} Pump;


typedef struct Chem_Data_structure
{
    int             NumVol;     /* Number of total volume in the rt simulator */
    int             NumOsv;     /* Number of grid blocks for the os3d (less ghost blocks) */
    int             NumFac;     /* Number of faces in the rt simulator        */
    int             NumDis;     /* Number of dispersive faces, used to calculate peclet number */
    int             NumEle;     /* Number of groundwater elements */
    int             NumRiv;     /* Number of river elements       */
    int             NumStc;     /* Number of total species in the rt simulator */
    int             NumSpc;     /* Number of primary species in the rt simulator      */
    int             NumSsc;     /* Number of secondary speices in the simulator */
    int             NumSdc;     /* Number of independent species (others depending on these species) */
    int             NumMin;     /* Number of minerals in the simulator */
    int             NumAds;     /* Number of adsorption species in the simulation */
    int             NumCex;     /* Number of cation exchange in the simulation */
    int             NumMkr;     /* Number of mineral kinetic reactions */
    int             NumAkr;     /* Number of aqueous kinetic reactions */
    int             OutItv;     /* controlling the output intervals, in unit of hours */
    int             TVDFlg;     /* TVD swith, 0 for off and 1 for on */
    int             SPCFlg;     /* speciation flg, 0 for total conc and 1 for pH */
    int             ACTmod;     /* activity coefficient mode, 0 for unity coefficient and 1 for DH equation */
    int             DHEdel;     /* update activity coefficient delay, 0 for delay ( update from previous time step) and 1 for solving the equation all together with speciation */
    int             TEMcpl;     /* whether or not to couple the temperature model */
    int             RelMin;     /* relative mineral flag. Determine whether the mineral volume fraction in input file is relative or absolution. For example, 0.1 calcite with RelMin =1 means the calcite occupy 10% of total solid volume. 0.1 calcite with RelMin = 0  means the calcite occupy 10% of the total pore volume. Be careful before using this key word */
    int             EffAds;     /* Keywords to control the usage of effective adsorption model */
    int             RecFlg;     /* A special flag to tell the code not do kinetic reaction, it will be much faster. Do not use it when you have kinetic reaction specified. Only suitable for doing stable isotope transport, etc. */
    int             PrpFlg;     /* Flag indicate whether or not the precipitation has been specified */
    int             AvgScl;     /* Flux will be averaged over a period of time, in the unit of min. E.g. AvgScl = 1 original flux from PIHM. AvgScl = 60, hourly flux updates from PIHM, it will be adaptively changed by the code itself.  */
    int             CptFlg;     /* Flag for coupling option. 1: coupling with pihm. 0: flow. others to be developped. */
    int             Delay;      /* RT start after PIHM running for a period of time, unit: days */
    int             React_delay;    /* Number of transport step per unit reaction step */
    int             NumPUMP;    /* Number of pumps  */
    int             SUFEFF;     /* surface effect */
    int             NumBTC;     /* Number of breakthrough points */
    int            *BTC_loc;    /* Array of locations of breakthrough points */
    int            *prepconcindex;  //
    double          TimLst;     /* Starting time step that is inherited from the last time step of hydro-reaction coupling */
    double          TimRiv;     /* transport time of river, calculated in fluxtrans and used in os3d. */
    double          Condensation;   /* a factor controls the concentration of infiltrating rain water as a ratio to the concentration in rain water */
    double          CnntVelo;   // velocity of minimum connected cells */
    double          CalPorosity;    // Porosity Calibration Coefficient, from Flux-PIHM */
    double          CalRate;    // 02.12 by Wei Zhi
    double          CalSSA;     // 02.12 by Wei Zhi
    double          CalGwinflux;    // 02.12 by Wei Zhi
    double          CalPrcpconc;    // 02.12 by Wei Zhi
    double          CalInitconc;    // 02.12 by Wei Zhi
    double          CalXsorption;   // 03.06 by Wei Zhi
    double        **Dependency; /* a matrix that describe the dependency of secondary species on the primary species */
    double        **Dep_kinetic;    /* a matrix that store the dependency of kinetic species on the primary species (sometimes conversion is required) */
    double        **Dep_kinetic_all;    // same as above, for all possible kinetic species
    double        **Totalconc;  /* a matrix that describe the contribution of each species to the total concentration */
#if NOT_YET_IMPLEMENTED
    double        **Totalconck; /* a matrix that describe the contribution of each species to the total concentration with kinetic reaction included */
#endif
    double         *Keq;        /* array of Keq s of the secondary species */
    double         *KeqKinect;  /* array of Keq s of the kinetic species */
    double         *KeqKinect_all;  /* same as above, all possible kinetic species */
    double          Temperature;    /* Temperature of the moment */
    double          StartTime;  /* start time of simulation, in unit of min   */
    double          Cementation;    /* Cementation factor, used to represent the connectivity of pores */
    double          rivd;       /* Stream discharge in cubic meter per day at the moment */
    double          riv;
    vol_conc       *Vcele;      // An array that stores the volumetric (vol) and chemical (conc) information of grid blocks
    vol_conc        Precipitation;  // The cell that stores the concentrations of chemicals in the rain.
    face           *Flux;       // connections between grid blocks
    species        *chemtype;   // information of chemical species
    Kinetic_Reaction *kinetics; // kinetics constants and dependencies.
    Debye_Huckel    DH;
    Pump           *pumps;      // injection/ groundwater contribution
    tsdata_struct  *TSD_prepconc;   // Time series data of concentration in precipitation.
} *Chem_Data;
