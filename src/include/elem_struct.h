#ifndef ELEM_STRUCT_HEADER
#define ELEM_STRUCT_HEADER

/* Element attribute */
typedef struct attrib_struct
{
    int             soil_type;              /* element soil type */
#ifdef _FBR_
    int             geol_type;              /* element geology type */
#endif
    int             lc_type;                /* element land cover type */
    int             bc_type[NUM_EDGE];      /* element boundary condition type
                                             */
#ifdef _FBR_
    int             fbrbc_type[NUM_EDGE];   /* element fractured bedrock layer
                                             * boundary condition type */
#endif
    int             meteo_type;             /* element meteorological forcing
                                             * type */
    int             lai_type;               /* element leaf area index forcing
                                             * type */
} attrib_struct;

/* Topography parameters */
typedef struct topo_struct
{
    double          area;                  /* area of element (m2) */
    double          x;                     /* x of centroid (m) */
    double          y;                     /* y of centroid (m) */
    double          zmin;                  /* soil bottom elevation (m) */
    double          zmax;                  /* surface elevation (m) */
#ifdef _FBR_
    double          zbed;                  /* impermeable bedrock elevation (m)
                                            */
#endif
    double          edge[NUM_EDGE];        /* length of edge (Edge i is from
                                            * node i to node i + 1) (m) */
    double          nabrdist[NUM_EDGE];    /* distance to neighbor (m) */
    double          nabr_x[NUM_EDGE];      /* x of neighbor centroid (m) */
    double          nabr_y[NUM_EDGE];      /* y of neighbor centroid (m) */
#ifdef _NOAH_
    double          slope;                 /* slope of element (degree) */
    double          aspect;                /* surface aspect of element (degree)
                                            */
    double          svf;                   /* sky view factor (-) */
    double          h_phi[36];             /* unobstructed angle in each
                                            * direction (degree) */
#endif
#ifdef _RT_
    double          areasub[NUM_EDGE];
#endif
} topo_struct;

/* Soil parameters */
typedef struct soil_struct
{
    double          depth;       /* soil depth (m) */
    double          ksath;       /* horizontal saturated hydraulic conductivity
                                  * (m s-1) */
    double          ksatv;       /* vertical saturated hydraulic conductivity
                                  * (m s-1) */
    double          kinfv;       /* saturated infiltration conductivity (m s-1)
                                  */
    double          dinf;        /* depth from ground surface across which head
                                  * gradient is calculated for infiltration (m)
                                  */
    double          alpha;       /* alpha from van Genuchten eqn (m-1) */
    double          beta;        /* beta (n) from van Genuchten eqn (-) */
    double          porosity;    /* soil porosity (m3 m-3) */
    double          smcmax;      /* maximum soil moisture content (m3 m-3) */
    double          smcmin;      /* residual soil moisture content (m3 m-3) */
    double          smcwlt;      /* wilting point (m3 m-3) */
    double          smcref;      /* soil moisture threshold where transpiration
                                  * begins to stress (m3 m-3) */
    double          dmac;        /* macropore depth (m) */
    double          kmach;       /* macropore horizontal saturated hydraulic
                                  * conductivity (m s-1) */
    double          kmacv;       /* macropore vertical saturated hydraulic
                                  * conductivity (m s-1) */
    double          areafv;      /* macropore area fraction on a vertical cross-
                                  * section (m2 m-2) */
    double          areafh;      /* macropore area fraction on a horizontal
                                  * cross-section (m2 m-2) */
#ifdef _NOAH_
    double          csoil;       /* soil heat capacity (J m-3 K-1) */
    double          quartz;      /* soil quartz content (-) */
    double          smcdry;      /* dry soil moisture threshold where direct
                                  * evap from top layer ends (m3 m-3) */
#endif
#ifdef _CYCLES_
    int             totalLayers;
    double          Curve_Number;
    double          Percent_Slope;
    double          annualTemperaturePhase;
    double          dampingDepth;
    double          cumulativeDepth[MAXLYR];
    double          nodeDepth[MAXLYR + 1];
    double          layerThickness[MAXLYR];
    double          Clay[MAXLYR];
    double          Sand[MAXLYR];
    double          IOM[MAXLYR];
    double          NO3[MAXLYR];
    double          NH4[MAXLYR];
    double          BD[MAXLYR];
    double          FC[MAXLYR];
    double          PWP[MAXLYR];
    double          Porosity[MAXLYR];
    double          PAW[MAXLYR];
    double          n2o[MAXLYR];
    double          SOC_Conc[MAXLYR];
    double          SOC_Mass[MAXLYR];
    double          SON_Mass[MAXLYR];
    double          MBC_Mass[MAXLYR];
    double          MBN_Mass[MAXLYR];
    double          SOCProfile;
    double          SONProfile;
    double          C_Humified;
    double          C_ResidueRespired;
    double          C_SoilRespired;
    double          soilTemperature[MAXLYR];
    double          waterContent[MAXLYR];
    double          waterUptake[MAXLYR];
    double          pH[MAXLYR];
    double          latflux[MAXLYR];
    double          evaporationVol;
    double          residueEvaporationVol;
    double          infiltrationVol;
    double          runoffVol;
    double          irrigationVol;
    double          drainageVol;
    double          NO3Leaching[NUM_EDGE];
    double          NH4Leaching[NUM_EDGE];
    double          NO3Profile;
    double          NH4Profile;
    double          N_Immobilization;
    double          N_Mineralization;
    double          N_NetMineralization;
    double          NH4_Nitrification;
    double          N2O_Nitrification;
    double          NO3_Denitrification;
    double          N2O_Denitrification;
    double          NH4_Volatilization;
    double          Rock[MAXLYR];
#endif
} soil_struct;

/* Fractured bedrock layer parameters */
typedef struct geol_struct
{
    double          depth;       /* soil depth (m) */
    double          ksath;       /* horizontal saturated hydraulic conductivity
                                  * (m s-1) */
    double          ksatv;       /* vertical saturated hydraulic conductivity
                                  * (m s-1) */
    double          alpha;       /* alpha from van Genuchten eqn (m-1) */
    double          beta;        /* beta (n) from van Genuchten eqn (-) */
    double          porosity;    /* porosity (m3 m-3) */
    double          smcmax;      /* maximum moisture content (m3 m-3) */
    double          smcmin;      /* residual moisture content (m3 m-3) */
} geol_struct;

/* Land cover parameters */
typedef struct lc_struct
{
    double          shdfac;       /* areal fractional coverage of green
                                   * vegetation (0.0-1.0) (-) */
    double          shdmin;       /* minimum areal fractional coverage of green
                                   * vegetation (-) */
    double          shdmax;       /* maximum areal fractional coverage of green
                                   * vegetation (-) */
    double          laimin;       /* minimum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    double          laimax;       /* maximum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    double          snup;         /* threshold snow depth (in water equivalent)
                                   * that implies 100% snow cover (m) */
    double          cfactr;       /* parameter used in the canopy interception
                                   * calculation (-) */
    double          emissmax;     /* minimum emissivity (-) */
    double          emissmin;     /* maximum emissivity (-) */
    double          albedomax;    /* minimum background albedo (-) */
    double          albedomin;    /* maximum background albedo (-) */
    double          z0max;        /* minimum roughness length (m) */
    double          z0min;        /* maximum roughness length (m) */
    double          rough;        /* surface roughness (Manning's n) (s m-1/3)
                                   */
    double          cmcfactr;     /* canopy water capacity per LAI (m) */
    int             bare;         /* flag that indicates bare ground */
    int             isurban;      /* flag that indicates urban */
} lc_struct;

/* Ecophysiological parameters */
typedef struct epconst_struct
{
    double          rsmin;                  /* minimum canopy resistance (s m-1)
                                             */
    double          rgl;                    /* reference incoming solar flux
                                             * for photosynthetically active
                                             * canopy (W m-2) */
    double          hs;                     /* parameter used in vapor pressure
                                             * deficit function (-) */
    double          topt;                   /* optimum transpiration air
                                             * temperature (K) */
    double          rsmax;                  /* cuticular resistance (s m-1) */
#ifdef _BGC_
    int             woody;                  /* flag: 1 = woody, 0 = non-woody */
    int             evergreen;              /* flag: 1 = evergreen,
                                             * 0 = deciduous */
    int             c3_flag;                /* flag: 1 = C3,  0 = C4 */
    int             phenology_flag;         /* flag: 1 = phenology mode
                                             * 0 = user defined */
    int             onday;                  /* day of year when leaves on */
    int             offday;                 /* day of year when leaves off */
    int             transfer_days;          /* growth period for transfer (day)
                                             */
    int             litfall_days;           /* growth period for litter fall
                                             * (day) */
    double          leaf_turnover;          /* annual leaf turnover fraction
                                             * (yr-1) */
    double          froot_turnover;         /* annual fine root turnover
                                             * fraction (yr-1) */
    double          livewood_turnover;      /* annual live wood turnover
                                             * fraction (yr-1) */
    double          daily_mortality_turnover;/* daily mortality turnover (day-1)
                                              */
    double          daily_fire_turnover;    /* daily fire turnover (day-1) */
    double          alloc_frootc_leafc;     /* new fine root C to new leaf C (-)
                                             */
    double          alloc_newstemc_newleafc;/* new stem C to new leaf C (-) */
    double          alloc_newlivewoodc_newwoodc;/* new livewood C:new wood C (-)
                                                 */
    double          alloc_crootc_stemc;     /* new live croot C to new live stem
                                             * C (-) */
    double          alloc_prop_curgrowth;   /* daily allocation to current
                                             * growth (-) */
    double          avg_proj_sla;           /* canopy average projected SLA
                                             * (m2 kgC-1) */
    double          sla_ratio;              /* ratio of shaded to sunlit
                                             * projected SLA (-) */
    double          lai_ratio;              /* ratio of (all-sided LA /
                                             * one-sided LA) (-) */
    double          ext_coef;               /* canopy light extinction
                                             * coefficient (-) */
    double          flnr;                   /* leaf N in Rubisco
                                             * (kgNRub kgNleaf-1) */
    double          psi_open;               /* psi at start of conductance
                                             * reduction (MPa) */
    double          psi_close;              /* psi at complete conductance
                                             * reduction (MPa) */
    double          vpd_open;               /* vpd at start of conductance
                                             * reduction (Pa) */
    double          vpd_close;              /* vpd at complete conductance
                                             * reduction (Pa) */
    double          froot_cn;               /* C:N for fine roots (kgC kgN-1) */
    double          leaf_cn;                /* C:N for leaves (kgC kgN-1) */
    double          livewood_cn;            /* C:N for live wood (kgC kgN-1) */
    double          deadwood_cn;            /* C:N for dead wood (kgC kgN-1) */
    double          leaflitr_cn;            /* constant C:N for leaf litter
                                             * (kgC kgN-1) */
    double          leaflitr_flab;          /* leaf litter labile fraction (-)
                                             */
    double          leaflitr_fucel;         /* leaf litter unshielded cellulose
                                             * fraction (-) */
    double          leaflitr_fscel;         /* leaf litter shielded cellulose
                                             * fraction (-) */
    double          leaflitr_flig;          /* leaf litter lignin fraction (-)
                                             */
    double          frootlitr_flab;         /* fine root litter labile fraction
                                             * (-) */
    double          frootlitr_fucel;        /* fine root litter unshielded
                                             * cellulose fraction (-) */
    double          frootlitr_fscel;        /* fine root litter shielded
                                             * cellulose fraction (-) */
    double          frootlitr_flig;         /* fine root litter lignin fraction
                                             * (-) */
    double          deadwood_fucel;         /* dead wood unshielded cellulose
                                             * fraction (-) */
    double          deadwood_fscel;         /* dead wood shielded cellulose
                                             * fraction (-) */
    double          deadwood_flig;          /* dead wood lignin fraction (-) */
#endif
} epconst_struct;

/* Physical states */
typedef struct pstate_struct
{
    double          rzd;                    /* rooting depth (m) */
    int             macpore_status;         /* macropore status */
    double          rc;                     /* canopy resistance (s m-1) */
    double          pc;                     /* plant coefficient (-) */
    double          proj_lai;               /* live projected leaf area index
                                             * (m2 m-2) */
    double          rcs;                    /* incoming solar rc factor (-) */
    double          rct;                    /* air temperature rc factor (-) */
    double          rcq;                    /* vapor pressure deficit rc factor
                                             * (-) */
    double          rcsoil;                 /* soil moisture rc factor (-) */
    double          albedo;                 /* surface albedo including snow
                                             * effect (-) */
    double          zlvl;                   /* height above ground of
                                             * atmospheric forcing variables (m)
                                             */
    double          zlvl_wind;              /* height above ground of wind
                                             * observations (m) */
    double          sfcspd;                 /* wind speed at height zlvl above
                                             * ground (m s-1) */
    double          rh;                     /* relative humidity (100%) */
    double          sfcprs;                 /* surface pressure at height zlvl
                                             * above ground (Pa) */
#ifdef _NOAH_
    double          alb;                    /* background snow-free surface
                                             * albedo (-) */
    double          snoalb;                 /* upper bound on maximum albedo
                                             * over deep snow (-) */
    int             nroot;                  /* number of root layers, a function
                                             * of vegetation type */
    double          rtdis[MAXLYR];          /* root distribution (-) */
    int             nsoil;                  /* number of soil layers */
    double          sldpth[MAXLYR];         /* thickness of each soil layer (m)
                                             */
    double          zsoil[MAXLYR];          /* distance from land surface to
                                             * bottom of each soil layer (m) */
    double          soilw;                  /* available soil moisture in root
                                             * zone (fraction between smcwlt and
                                             * smcmax) (-) */
    double          frzk;                   /* frozen ground parameter (-) */
    double          frzx;                   /* adjusted frozen ground parameter
                                             * (-) */
    double          czil;                   /* Zilitinkevich constant (-) */
    double          emissi;                 /* surface emissivity (between 0 and
                                             * 1) (-) */
    double          ch;                     /* surface exchange coefficient for
                                             * heat and moisture (m s-1) */
    double          cm;                     /* surface exchange coefficient for
                                             * momentum (m s-1) */
    double          rch;                    /* = ch * air density * CP
                                             * (W m-2 K-1) */
    double          z0;                     /* time varying roughness length as
                                             * function of snow depth (-) */
    double          fcr;                    /* reduction of infiltration caused
                                             * by frozen ground (-) */
    int             nmacd;                  /* number of soil layers with
                                             * macropore */
    double          salp;                   /* shape parameter of distribution
                                             * function of snow cover (-) */
    double          fxexp;                  /* soil evaporation exponent used
                                             * in direct evaporation (-) */
    double          sbeta;                  /* parameter used to calculate
                                             * vegetation effect on soil heat
                                             * (-) */
    double          lvcoef;                 /* parameter controls surface snow
                                             * albedo in the presence of snow
                                             * cover (-) */
    double          snotime1;               /* age of the snow on the ground (s)
                                             */
    double          ribb;                   /* bulk Richardson number used to
                                             * limit the dew/frost (-) */
    double          beta;                   /* ratio of actual/potential evap
                                             * (-) */
    double          sncovr;                 /* fractional snow cover (-) */
    double          q1;                     /* effective mixing ratio at surface
                                             * (kg kg-1) */
    double          q2;                     /* mixing ratio at height zlvl above
                                             * (kg kg-1) */
    double          ffrozp;                 /* fraction of frozen precipitation
                                             * (-) */
    double          z0brd;                  /* background fixed roughness length
                                             * (-) */
    double          embrd;                  /* background surface emissivity
                                             * (-) */
    double          q2sat;                  /* saturation air humidity at height
                                             * zlvl above ground (kg kg-1) */
    double          q2d;                    /* air humidity deficit (kg kg-1) */
    double          dqsdt2;                 /* slope of saturation specific
                                             * humidity curve at T = sfctmp
                                             * (kg kg-1 K-1) */
    int             nwtbl;                  /* layer where water table is within
                                             */
    double          sndens;                 /* snow density (dimensionless
                                             * fraction of H2O density) (-) */
    double          snowh;                  /* actual snow depth (-) */
    double          sncond;                 /* snow thermal conductivity
                                             * (W m-1 K-1) */
    double          rr;                     /* parameter in Penman potential
                                             * evaporation (-) */
    double          epsca;                  /* parameter in Penman potential
                                             * evaporation (K) */
    double          eta_kinematic;          /* actual latent heat flux
                                             * (kg m-2 s-1) */
    double          zbot;                   /* depth of lower boundary soil
                                             * temperature (m) */
    double          tbot;                   /* bottom soil temperature (local
                                             * yearly-mean sfc air temperature)
                                             * (K) */
    double          gwet;                   /* fraction of transpiration from
                                             * groundwater (-) */
    double          satdpth[MAXLYR];        /* depth of groundwater in each soil
                                             * layer (m) */
#endif
#ifdef _BGC_
    double          co2;                    /* atmospheric CO2 concentration
                                             * (ppm) */
    double          ppfd_per_plaisun;       /* ppfd per unit sunlit proj LAI
                                             * (umol m-2 s-1) */
    double          ppfd_per_plaishade;     /* ppfd per unit shaded proj LAI
                                             * (umol m-2 s-1) */
    double          all_lai;                /* live all-sided leaf area index
                                             * (m2 m-2) */
    double          plaisun;                /* sunlit projected leaf area index
                                             * (m2 m-2) */
    double          plaishade;              /* shaded projected leaf area index
                                             * (m2 m-2) */
#endif
} pstate_struct;

/* Water states */
typedef struct wstate_struct
{
    double          surf;            /* equivalent surface water level (m) */
    double          unsat;           /* unsaturated zone water storage (m) */
    double          gw;              /* groundwater level (m) */
    double          sneqv;           /* liquid water-equivalent snow depth (m)
                                      */
    double          cmcmax;          /* maximum canopy water capacity (m) */
    double          cmc;             /* interception storage (m) */
    double          surfh;           /* actual surface water level (m) */
#ifdef _FBR_
    double          fbr_unsat;       /* unsaturated storage in fractured bedrock
                                      * layer (m) */
    double          fbr_gw;          /* deep groundwater in fractured bedrock
                                      * layer (m) */
#endif
#ifdef _NOAH_
    double          smc[MAXLYR];     /* total soil moisture content (m3 m-3) */
    double          sh2o[MAXLYR];    /* unfrozen soil moisture content (m3 m-3)
                                      */
    double          soilm;           /* total soil column moisture content (m)
                                      */
#endif
} wstate_struct;

/* Water fluxes */
typedef struct wflux_struct
{
    double          ovlflow[NUM_EDGE];      /* overland flow (m3 s-1) */
    double          subsurf[NUM_EDGE];      /* subsurface flow (m3 s-1) */
    double          prcp;                   /* precipitation on each element
                                             * (m s-1) */
    double          pcpdrp;                 /* combined prcp and drip (from
                                             * canopy) that goes into the soil
                                             * (m s-1) */
    double          infil;                  /* variable infiltration rate
                                             * (m s-1) */
    double          rechg;                  /* recharge rate to groundwater
                                             * (m s-1) */
    double          drip;                   /* through-fall of precipitation
                                             * and/or dew (m s-1) */
    double          edir;                   /* direct soil evaporation (m s-1)
                                             */
    double          ett;                    /* total plant transpiration (m s-1)
                                             */
    double          ec;                     /* canopy water evaporation (m s-1)
                                             */
    double          etp;                    /* potential evaporation (m s-1) */
    double          eta;                    /* actual evapotranspiration (m s-1)
                                             */
    double          edir_surf;              /* direct evaporation from surface
                                             * water (m s-1) */
    double          edir_unsat;             /* direct evaporation from
                                             * unsaturated zone (m s-1) */
    double          edir_gw;                /* direct evaporation from saturated
                                             * zone (m s-1) */
    double          ett_unsat;              /* transpiration from unsaturated
                                             * zone (m s-1) */
    double          ett_gw;                 /* transpiration from saturated zone
                                             * (m s-1) */
    double          esnow;                  /* sublimation from (or deposition
                                             * to) snowpack (m s-1); */
#ifdef _FBR_
    double          fbr_infil;              /* fractured bedrock infiltration
                                             * (m s-1) */
    double          fbr_rechg;              /* fractured bedrock recharge
                                             * (m s-1) */
    double          fbrflow[NUM_EDGE];      /* lateral fractured bedrock flow
                                             * (m3 s-1) */
#endif
#ifdef _NOAH_
    double          et[MAXLYR];             /* plant transpiration from each
                                             * soil layer (m s-1) */
    double          runoff2;                /* total subsurface flow (m s-1) */
    double          runoff2_lyr[MAXLYR];    /* subsurface flow from each soil
                                             * layer (m s-1) */
    double          runoff3;                /* numerical truncation in excess
                                             * of porosity (smcmax) for a given
                                             * soil layer at the end of a time
                                             * step (m s-1) */
    double          smflxv[MAXLYR];         /* vertical soil moisture flux
                                             * between soil layers (m s-1) */
    double          smflxh[NUM_EDGE][MAXLYR];/* horizontal soil moisture flux at
                                             * each soil layer from each edge
                                             * (m s-1) */
    double          dew;                    /* dewfall (or frostfall for
                                             * T < 273.15) (m s-1) */
    double          snomlt;                 /* water equivalent snow melt
                                             * (m s-1) */
    double          etns;                   /* (m s-1) */
#endif
#ifdef _CYCLES_
    double          eres;                   /* evaporation from residue (m s-1)
                                             */
#endif
} wflux_struct;

/* Energy states */
typedef struct estate_struct
{
    double          sfctmp;         /* air temperature at height zlvl above
                                     * ground (K) */
#ifdef _NOAH_
    double          t1;             /* ground/canopy/snowpack effective skin
                                     * temperature (K) */
    double          th2;            /* air potential temperature at height zlvl
                                     * above ground (K) */
    double          stc[MAXLYR];    /* soil temperature (K) */
#endif
} estate_struct;

/* Energy fluxes */
typedef struct eflux_struct
{
    double          soldn;                  /* solar downward radiation (W m-2)
                                             */
#ifdef _NOAH_
    double          solnet;                 /* net downward solar radiation
                                             * (W m-2) */
    double          etp;                    /* potential evaporation (W m-2) */
    double          ssoil;                  /* soil heat flux (W m-2) */
    double          eta;                    /* actual latent heat flux (W m-2)
                                             */
    double          sheat;                  /* sensible heat flux (W m-2) */
    double          fdown;                  /* radiation forcing at the surface
                                             * (W m-2) */
    double          lwdn;                   /* absorbed longwave downward
                                             * radiation (W m-2) */
    double          ec;                     /* canopy water evaporation (W m-2)
                                             */
    double          edir;                   /* direct soil evaporation (W m-2)
                                             */
    double          et[MAXLYR];             /* plant transpiration from each
                                             * soil layer (W m-2) */
    double          ett;                    /* total plant transpiration (W m-2)
                                             */
    double          esnow;                  /* sublimation from (or deposition)
                                             * snowpack (W m-2) */
    double          soldir;                 /* direct solar radiation (W m-2) */
    double          soldif;                 /* diffused solar radiation (W m-2)
                                             */
    double          longwave;               /* longwave radiation forcing
                                             * (W m-2) */
    double          flx1;                   /* latent heat flux from
                                             * precipitation accumulating as
                                             * snow (W m-2) */
    double          flx2;                   /* freezing rain latent heat flux
                                             * (W m-2) */
    double          flx3;                   /* snow melt latent heat flux
                                             * (W m-2) */
#endif
#ifdef _BGC_
    double          swabs_per_plaisun;      /* swabs per unit sunlit proj LAI
                                             * (W m-2) */
    double          swabs_per_plaishade;    /* swabs per unit shaded proj LAI
                                             * (W m-2) */
#endif
} eflux_struct;

#ifdef _CYCLES_
typedef struct crop_struct
{
    /* Instance of the crop that is being planted */
    char            cropName[MAXSTRING];

    /* User Defined Auto Irrigation */
    int             autoIrrigationUsed;
    int             autoIrrigationStartDay;
    int             autoIrrigationStopDay;
    double          autoIrrigationWaterDepletion;
    int             autoIrrigationLastSoilLayer;

    /* User Defined Auto Fertilization */
    int             autoFertilizationUsed;
    int             autoFertilizationStartDay;
    int             autoFertilizationStopDay;
    double          autoFertilizationMass;
    char            autoFertilizationSource;
    char            autoFertilizationForm;
    int             autoFertilizationMethod;

    /* Crop Status Flags */
    int             cropGrowing;
    int             cropMature;

    /* State Variables */
    double          svTT_Daily;
    double          svTT_Cumulative;
    double          svRadiationInterception;
    double          svBiomass;
    double          svShoot;
    double          svRoot;
    double          svRizho;
    double          svShootDailyGrowth;
    double          svRootDailyGrowth;
    double          svRizhoDailyDeposition;
    double          svUnstressedShootDailyGrowth;
    double          svUnstressedRootDailyGrowth;
    double          svPostFloweringShootBiomass;
    double          svRootingDepth;
    double          svTranspiration;
    double          svTranspirationPotential;
    double          svN_Shoot;
    double          svN_Root;
    double          svN_Rhizo;
    double          svN_RizhoDailyDeposition;
    double          svN_AutoAdded;
    double          svN_Fixation;
    double          svWaterStressFactor;
    double          svN_StressFactor;

    double          svShootUnstressed;
    double          svN_StressCumulative;

    double          svRadiationInterception_nc;

    double          dailyTranspiration;
    double          dailyTranspirationPotential;

    double          userFloweringTT;
    double          userMaturityTT;
    double          userMaximumSoilCoverage;
    double          userMaximumRootingDepth;
    double          userExpectedYieldAvg;
    double          userExpectedYieldMax;
    double          userExpectedYieldMin;
    double          userPercentMoistureInYield;
    double          userFractionResidueStanding;
    double          userFractionResidueRemoved;
    double          userClippingBiomassThresholdUpper;
    double          userClippingBiomassThresholdLower;
    double          userClippingTiming;
    int             userClippingDestiny;
    double          userTranspirationMinTemperature;
    double          userTranspirationThresholdTemperature;
    double          userColdDamageMinTemperature;
    double          userColdDamageThresholdTemperature;
    double          userTemperatureBase;
    double          userTemperatureOptimum;
    double          userTemperatureMaximum;
    double          userShootPartitionInitial;
    double          userShootPartitionFinal;
    double          userRadiationUseEfficiency;
    double          userTranspirationUseEfficiency;
    double          userHIx;
    double          userHIo;    /* intercept harvest index */
    double          userHIk;
    double          userEmergenceTT;
    double          userNMaxConcentration;
    double          userNDilutionSlope;
    double          userKc;
    int             userAnnual;
    int             userLegume;
    int             userC3orC4;
    double          userExtinctionCoefficient;

    double          userPlantingDensity;

    int             userClippingStart;
    int             userClippingEnd;
    double          calculatedSimAvgYield;
    double          calculatedSimMaxYield;
    double          calculatedSimMinYield;
    double          LWP_StressOnset;
    double          LWP_WiltingPoint;
    double          transpirationMax;

    int             harvestDateFinal;
    int             harvestCount;
    int             stageGrowth;

    double          rcForageYield;
    double          rcGrainYield;
    double          rcBiomass;
    double          rcRoot;
    double          rcResidueBiomass;
    double          rcCropTranspiration;
    double          rcCropTranspirationPotential;
    double          rcSoilWaterEvaporation;
    double          rcHarvestIndex;
    double          rcTotalNitrogen;
    double          rcRootNitrogen;
    double          rcGrainNitrogenYield;
    double          rcForageNitrogenYield;
    double          rcNitrogenCumulative;
    double          rcNitrogenInHarvest;
    double          rcNitrogenInResidue;
    double          rcNitrogenForageConc;
} crop_struct;

typedef struct comm_struct
{
    /* State Variables */
    double          svRadiationInterception;
    double          svBiomass;
    double          svShoot;
    double          svRoot;
    double          svRizho;
    double          svShootDailyGrowth;
    double          svRootDailyGrowth;
    double          svRizhoDailyDeposition;
    double          svRootingDepth; /* maximum */
    double          svTranspiration;
    double          svTranspirationPotential;
    double          svN_Shoot;
    double          svN_Root;
    double          svN_Rhizo;
    double          svN_RizhoDailyDeposition;
    double          svN_AutoAdded;
    double          svN_Fixation;
    double          svWaterStressFactor;
    double          svN_StressFactor;

    crop_struct    *Crop;
    int             NumCrop;
    int             NumActiveCrop;
} comm_struct;

typedef struct residue_struct
{
    double          residueInterception;
    double          stanResidueTau;
    double          flatResidueTau;
    double          stanResidueMass;
    double          flatResidueMass;
    double          stanResidueN;
    double          flatResidueN;
    double          manureSurfaceC;
    double          manureSurfaceN;

    double          stanResidueWater;
    double          flatResidueWater;   /* (mm) */

    double          residueAbgd[MAXLYR];
    double          residueRt[MAXLYR];
    double          residueRz[MAXLYR];
    double          residueAbgdN[MAXLYR];
    double          residueRtN[MAXLYR];
    double          residueRzN[MAXLYR];
    double          yearResidueBiomass;
    double          yearResidueHarvested;
    double          yearRootBiomass;
    double          yearRhizodepositionBiomass;
    double          manureC[MAXLYR];
    double          manureN[MAXLYR];    /* Mg/ha */
} residue_struct;

typedef struct soilc_struct
{
    double          factorComposite[MAXLYR];
    double          carbonRespired[MAXLYR];
    double          rootBiomassInput[MAXLYR];
    double          rhizBiomassInput[MAXLYR];
    double          abgdBiomassInput[MAXLYR];
    double          rootCarbonInput[MAXLYR];
    double          rhizCarbonInput[MAXLYR];
    double          manuCarbonInput[MAXLYR];
    double          abgdCarbonInput[MAXLYR];
    double          carbonMassInitial[MAXLYR];
    double          carbonMassFinal[MAXLYR];
    double          annualDecompositionFactor[MAXLYR];
    double          annualSoilCarbonDecompositionRate[MAXLYR];
    double          annualCarbonInputByLayer[MAXLYR];
    double          annualHumifiedCarbonMass[MAXLYR];
    double          annualRespiredCarbonMass[MAXLYR];
    double          annualRespiredResidueCarbonMass[MAXLYR];
    double          annualHumificationCoefficient[MAXLYR];
    double          annualNmineralization[MAXLYR];
    double          annualNImmobilization[MAXLYR];
    double          annualNNetMineralization[MAXLYR];
    double          annualAmmoniumNitrification;
    double          annualNitrousOxidefromNitrification;
    double          annualAmmoniaVolatilization;
    double          annualNO3Denitrification;
    double          annualNitrousOxidefromDenitrification;
    double          annualNitrateLeaching;
    double          annualAmmoniumLeaching;
} soilc_struct;

typedef struct cyclesic_struct
{
    double          SOC_Mass[MAXLYR];
    double          SON_Mass[MAXLYR];
    double          MBC_Mass[MAXLYR];
    double          MBN_Mass[MAXLYR];
    double          NO3[MAXLYR];
    double          NH4[MAXLYR];
} cyclesic_struct;
#endif

/* Boundary conditions */
typedef struct bc_struct
{
    double          head[NUM_EDGE];    /* value of Dirichlet-type boundary
                                        * condition (m) */
    double          flux[NUM_EDGE];    /* value of Neumann-type boundary
                                        * condition (m3 s-1) */
} bc_struct;

/* Land surface and hydrologic initial conditions */
typedef struct ic_struct
{
    double          cmc;
    double          sneqv;
    double          surf;
    double          unsat;
    double          gw;
#ifdef _FBR_
    double          fbr_unsat;
    double          fbr_gw;
#endif
#ifdef _NOAH_
    double          t1;
    double          snowh;
    double          stc[MAXLYR];
    double          smc[MAXLYR];
    double          sh2o[MAXLYR];
#endif
} ic_struct;

#ifdef _BGC_
/* CN initial conditions */
typedef struct bgcic_struct
{
    double          leafc;
    double          leafc_storage;
    double          leafc_transfer;
    double          frootc;
    double          frootc_storage;
    double          frootc_transfer;
    double          livestemc;
    double          livestemc_storage;
    double          livestemc_transfer;
    double          deadstemc;
    double          deadstemc_storage;
    double          deadstemc_transfer;
    double          livecrootc;
    double          livecrootc_storage;
    double          livecrootc_transfer;
    double          deadcrootc;
    double          deadcrootc_storage;
    double          deadcrootc_transfer;
    double          gresp_storage;
    double          gresp_transfer;
    double          cwdc;
    double          litr1c;
    double          litr2c;
    double          litr3c;
    double          litr4c;
    double          soil1c;
    double          soil2c;
    double          soil3c;
    double          soil4c;
    double          cpool;
    double          leafn;
    double          leafn_storage;
    double          leafn_transfer;
    double          frootn;
    double          frootn_storage;
    double          frootn_transfer;
    double          livestemn;
    double          livestemn_storage;
    double          livestemn_transfer;
    double          deadstemn;
    double          deadstemn_storage;
    double          deadstemn_transfer;
    double          livecrootn;
    double          livecrootn_storage;
    double          livecrootn_transfer;
    double          deadcrootn;
    double          deadcrootn_storage;
    double          deadcrootn_transfer;
    double          cwdn;
    double          litr1n;
    double          litr2n;
    double          litr3n;
    double          litr4n;
    double          soil1n;
    double          soil2n;
    double          soil3n;
    double          soil4n;
    double          surfn;
    double          sminn;
    double          retransn;
    double          npool;
    double          prev_leafc_to_litter;
    double          prev_frootc_to_litter;
    double          dsr;
    int             dormant_flag;
    int             onset_flag;
    int             onset_counter;
    int             onset_gddflag;
    double          onset_fdd;
    double          onset_gdd;
    double          onset_swi;
    int             offset_flag;
    int             offset_counter;
    double          offset_fdd;
    double          offset_swi;
} bgcic_struct;

/* CN spinup variables */
typedef struct spinup_struct
{
    double          soilc_prev;
    double          totalc_prev;
    double          soilc;
    double          totalc;
    int             steady;
} spinup_struct;
#endif

#ifdef _CYCLES_
typedef struct weather_struct
{
    int            *lastDoy;
    double        **wind;
    double        **solarRadiation;
    double        **tMax;
    double        **tMin;
    double        **vpd;
    double          atmosphericPressure;
} weather_struct;

typedef struct snow_struct
{
    double          snowCover;
} snow_struct;

typedef struct solute_struct
{
    double          soluteMass[MAXLYR];
    double          soluteMassAdsorbed[MAXLYR];
    double          soluteConc[MAXLYR];
    double          soluteFluxLat[NUM_EDGE][MAXLYR];
    double          soluteFluxVert[MAXLYR];
} solute_struct;
#endif

#ifdef _DAILY_

/* Daily average variables */
typedef struct daily_struct
{
    int             counter;               /* counter used for averaging */
    int             daylight_counter;      /* counter used for daytime averaging
                                            */
    double          avg_surf;              /* daily average surface water level
                                            * (m) */
    double          avg_unsat;             /* daily average unsaturated water
                                            * storage (m) */
    double          avg_gw;                /* daily average groundwater level
                                            * (m) */
    double          avg_sh2o[MAXLYR];      /* daily average unfrozen soil water
                                            * content (m3 m-3) */
    double          avg_smc[MAXLYR];       /* daily average unfrozen soil
                                            * moisture content (m3 m-3) */
    double          avg_smflxv[MAXLYR];    /* daily average vertical soil
                                            * moisture flux (m s-1) */
    double          avg_et[MAXLYR];        /* daily average evapotranspiration
                                            * (m s-1) */
    double          avg_q2d;               /* daily average mixing ratio deficit
                                            * (kg kg-1) */
    double          avg_sfcprs;            /* daily average air pressure (Pa) */
    double          avg_ch;                /* daily average surface exchange
                                            * coefficient (m s-1) */
    double          avg_rc;                /* daily average stomatal resistance
                                            * (s m-1) */
    double          avg_albedo;            /* daily average surface albedo
                                            * (-) */
    double          avg_sfcspd;            /* daily average wind speed (m s-1)
                                            */
    double          avg_sncovr;            /* daily average snow cover fraction
                                            * (-) */
    double          tmax;                  /* daily maximum air temperature (K)
                                            */
    double          tmin;                  /* daily minimum air temperature (K)
                                            */
    double          avg_sfctmp;            /* daily average air temperature (K)
                                            */
    double          tday;                  /* daytime average air temperature
                                            * (K) */
    double          tnight;                /* nighttime average air temperature
                                            * (K) */
    double          avg_stc[MAXLYR];       /* daily average soil temperature (K)
                                            */
    double          avg_soldn;             /* daytime average downward solar
                                            * radiation (W m-2) */
    double          solar_total;           /* daily total solar radiation
                                            * (J m-2) */
} daily_struct;
#endif

#ifdef _BGC_
/* Carbon state variables (including sums for sources and sinks) */
typedef struct cstate_struct
{
    double          leafc;                  /* leaf C (kgC m-2) */
    double          leafc_storage;          /* leaf C storage (kgC m-2) */
    double          leafc_transfer;         /* leaf C transfer (kgC m-2) */
    double          frootc;                 /* fine root C (kgC m-2) */
    double          frootc_storage;         /* fine root C storage (kgC m-2) */
    double          frootc_transfer;        /* fine root C transfer (kgC m-2) */
    double          livestemc;              /* live stem C (kgC m-2) */
    double          livestemc_storage;      /* live stem C storage (kgC m-2) */
    double          livestemc_transfer;     /* live stem C transfer (kgC m-2) */
    double          deadstemc;              /* dead stem C (kgC m-2) */
    double          deadstemc_storage;      /* dead stem C storage (kgC m-2) */
    double          deadstemc_transfer;     /* dead stem C transfer (kgC m-2) */
    double          livecrootc;             /* live coarse root C (kgC m-2) */
    double          livecrootc_storage;     /* live coarse root C storage
                                             * (kgC m-2) */
    double          livecrootc_transfer;    /* live coarse root C transfer
                                             * (kgC m-2) */
    double          deadcrootc;             /* dead coarse root C (kgC m-2) */
    double          deadcrootc_storage;     /* dead coarse root C storage
                                             * (kgC m-2) */
    double          deadcrootc_transfer;    /* dead coarse root C transfer
                                             * (kgC m-2) */
    double          gresp_storage;          /* growth respiration storage
                                             * (kgC m-2) */
    double          gresp_transfer;         /* growth respiration transfer
                                             * (kgC m-2) */
    double          cwdc;                   /* coarse woody debris C (kgC m-2)
                                             */
    double          litr1c;                 /* litter labile C (kgC m-2) */
    double          litr2c;                 /* litter unshielded cellulose C
                                             * (kgC m-2) */
    double          litr3c;                 /* litter shielded cellulose C
                                             * (kgC m-2) */
    double          litr4c;                 /* litter lignin C (kgC m-2) */
    double          soil1c;                 /* microbial recycling pool C (fast)
                                             * (kgC m-2) */
    double          soil2c;                 /* microbial recycling pool C
                                             * (medium)
                                             * (kgC m-2) */
    double          soil3c;                 /* microbial recycling pool C (slow)
                                             * (kgC m-2) */
    double          soil4c;                 /* recalcitrant SOM C (humus,
                                             * slowest) (kgC m-2) */
    double          cpool;                  /* temporary photosynthate C pool
                                             * (kgC m-2) */
    double          psnsun_src;             /* SUM of gross PSN from sunlit
                                             * canopy (kgC m-2) */
    double          psnshade_src;           /* SUM of gross PSN from shaded
                                             * canopy (kgC m-2) */
    double          leaf_mr_snk;            /* SUM of leaf maint resp (kgC m-2)
                                             */
    double          leaf_gr_snk;            /* SUM of leaf growth resp (kgC m-2)
                                             */
    double          froot_mr_snk;           /* SUM of fine root maint resp
                                             * (kgC m-2) */
    double          froot_gr_snk;           /* SUM of fine root growth resp
                                             * (kgC m-2) */
    double          livestem_mr_snk;        /* SUM of live stem maint resp
                                             * (kgC m-2) */
    double          livestem_gr_snk;        /* SUM of live stem growth resp
                                             * (kgC m-2) */
    double          deadstem_gr_snk;        /* SUM of dead stem growth resp
                                             * (kgC m-2) */
    double          livecroot_mr_snk;       /* SUM of live coarse root maint
                                             * resp (kgC m-2) */
    double          livecroot_gr_snk;       /* SUM of live coarse root growth
                                             * resp (kgC m-2) */
    double          deadcroot_gr_snk;       /* SUM of dead coarse root growth
                                             * resp (kgC m-2) */
    double          litr1_hr_snk;           /* SUM of labile litr microbial resp
                                             * (kgC m-2) */
    double          litr2_hr_snk;           /* SUM of cellulose litr microbial
                                             * resp (kgC m-2) */
    double          litr4_hr_snk;           /* SUM of lignin litr microbial resp
                                             * (kgC m-2) */
    double          soil1_hr_snk;           /* SUM of fast microbial respiration
                                             * (kgC m-2) */
    double          soil2_hr_snk;           /* SUM of medium microbial
                                             * respiration (kgC m-2) */
    double          soil3_hr_snk;           /* SUM of slow microbial respiration
                                             * (kgC m-2) */
    double          soil4_hr_snk;           /* SUM of recalcitrant SOM
                                             * respiration (kgC m-2) */
    double          fire_snk;               /* SUM of fire losses (kgC m-2) */
} cstate_struct;

/* Daily carbon flux variables */
typedef struct cflux_struct
{
    /* mortality fluxes (kgC m-2 day-1) */
    double          m_leafc_to_litr1c;
    double          m_leafc_to_litr2c;
    double          m_leafc_to_litr3c;
    double          m_leafc_to_litr4c;
    double          m_frootc_to_litr1c;
    double          m_frootc_to_litr2c;
    double          m_frootc_to_litr3c;
    double          m_frootc_to_litr4c;
    double          m_leafc_storage_to_litr1c;
    double          m_frootc_storage_to_litr1c;
    double          m_livestemc_storage_to_litr1c;
    double          m_deadstemc_storage_to_litr1c;
    double          m_livecrootc_storage_to_litr1c;
    double          m_deadcrootc_storage_to_litr1c;
    double          m_leafc_transfer_to_litr1c;
    double          m_frootc_transfer_to_litr1c;
    double          m_livestemc_transfer_to_litr1c;
    double          m_deadstemc_transfer_to_litr1c;
    double          m_livecrootc_transfer_to_litr1c;
    double          m_deadcrootc_transfer_to_litr1c;
    double          m_livestemc_to_cwdc;
    double          m_deadstemc_to_cwdc;
    double          m_livecrootc_to_cwdc;
    double          m_deadcrootc_to_cwdc;
    double          m_gresp_storage_to_litr1c;
    double          m_gresp_transfer_to_litr1c;
    /* fire fluxes (kgC m-2 day-1) */
    double          m_leafc_to_fire;
    double          m_frootc_to_fire;
    double          m_leafc_storage_to_fire;
    double          m_frootc_storage_to_fire;
    double          m_livestemc_storage_to_fire;
    double          m_deadstemc_storage_to_fire;
    double          m_livecrootc_storage_to_fire;
    double          m_deadcrootc_storage_to_fire;
    double          m_leafc_transfer_to_fire;
    double          m_frootc_transfer_to_fire;
    double          m_livestemc_transfer_to_fire;
    double          m_deadstemc_transfer_to_fire;
    double          m_livecrootc_transfer_to_fire;
    double          m_deadcrootc_transfer_to_fire;
    double          m_livestemc_to_fire;
    double          m_deadstemc_to_fire;
    double          m_livecrootc_to_fire;
    double          m_deadcrootc_to_fire;
    double          m_gresp_storage_to_fire;
    double          m_gresp_transfer_to_fire;
    double          m_litr1c_to_fire;
    double          m_litr2c_to_fire;
    double          m_litr3c_to_fire;
    double          m_litr4c_to_fire;
    double          m_cwdc_to_fire;
    /* phenology fluxes from transfer pools (kgC m-2 day-1) */
    double          leafc_transfer_to_leafc;
    double          frootc_transfer_to_frootc;
    double          livestemc_transfer_to_livestemc;
    double          deadstemc_transfer_to_deadstemc;
    double          livecrootc_transfer_to_livecrootc;
    double          deadcrootc_transfer_to_deadcrootc;
    /* leaf and fine root litterfall (kgC m-2 day-1) */
    double          leafc_to_litr1c;
    double          leafc_to_litr2c;
    double          leafc_to_litr3c;
    double          leafc_to_litr4c;
    double          frootc_to_litr1c;
    double          frootc_to_litr2c;
    double          frootc_to_litr3c;
    double          frootc_to_litr4c;
    /* maintenance respiration fluxes (kgC m-2 day-1) */
    double          leaf_day_mr;
    double          leaf_night_mr;
    double          froot_mr;
    double          livestem_mr;
    double          livecroot_mr;
    /* photosynthesis fluxes (kgC m-2 day-1) */
    double          psnsun_to_cpool;
    double          psnshade_to_cpool;
    /* litter decomposition fluxes (kgC m-2 day-1) */
    double          cwdc_to_litr2c;
    double          cwdc_to_litr3c;
    double          cwdc_to_litr4c;
    double          litr1_hr;
    double          litr1c_to_soil1c;
    double          litr2_hr;
    double          litr2c_to_soil2c;
    double          litr3c_to_litr2c;
    double          litr4_hr;
    double          litr4c_to_soil3c;
    double          soil1_hr;
    double          soil1c_to_soil2c;
    double          soil2_hr;
    double          soil2c_to_soil3c;
    double          soil3_hr;
    double          soil3c_to_soil4c;
    double          soil4_hr;
    /* daily allocation fluxes from current GPP (kgC m-2 day-1) */
    double          cpool_to_leafc;
    double          cpool_to_leafc_storage;
    double          cpool_to_frootc;
    double          cpool_to_frootc_storage;
    double          cpool_to_livestemc;
    double          cpool_to_livestemc_storage;
    double          cpool_to_deadstemc;
    double          cpool_to_deadstemc_storage;
    double          cpool_to_livecrootc;
    double          cpool_to_livecrootc_storage;
    double          cpool_to_deadcrootc;
    double          cpool_to_deadcrootc_storage;
    double          cpool_to_gresp_storage;
    /* daily growth respiration fluxes (kgC m-2 day-1) */
    double          cpool_leaf_gr;
    double          cpool_leaf_storage_gr;
    double          transfer_leaf_gr;
    double          cpool_froot_gr;
    double          cpool_froot_storage_gr;
    double          transfer_froot_gr;
    double          cpool_livestem_gr;
    double          cpool_livestem_storage_gr;
    double          transfer_livestem_gr;
    double          cpool_deadstem_gr;
    double          cpool_deadstem_storage_gr;
    double          transfer_deadstem_gr;
    double          cpool_livecroot_gr;
    double          cpool_livecroot_storage_gr;
    double          transfer_livecroot_gr;
    double          cpool_deadcroot_gr;
    double          cpool_deadcroot_storage_gr;
    double          transfer_deadcroot_gr;
    /* annual turnover of storage to transfer pools (kgC m-2 day-1) */
    double          leafc_storage_to_leafc_transfer;
    double          frootc_storage_to_frootc_transfer;
    double          livestemc_storage_to_livestemc_transfer;
    double          deadstemc_storage_to_deadstemc_transfer;
    double          livecrootc_storage_to_livecrootc_transfer;
    double          deadcrootc_storage_to_deadcrootc_transfer;
    double          gresp_storage_to_gresp_transfer;
    /* turnover of live wood to dead wood (kgC m-2 day-1) */
    double          livestemc_to_deadstemc;
    double          livecrootc_to_deadcrootc;
} cflux_struct;


/* Nitrogen state variables (including sums for sources and sinks) */
typedef struct nstate_struct
{
    /* leaf N (kgN m-2) */
    double          leafn;
    double          leafn_storage;
    double          leafn_transfer;
    /* fine root N (kgN m-2) */
    double          frootn;
    double          frootn_storage;
    double          frootn_transfer;
    /* live stem N (kgN m-2) */
    double          livestemn;
    double          livestemn_storage;
    double          livestemn_transfer;
    /* dead stem N (kgN m-2) */
    double          deadstemn;
    double          deadstemn_storage;
    double          deadstemn_transfer;
    /* live coarse root N (kgN m-2) */
    double          livecrootn;
    double          livecrootn_storage;
    double          livecrootn_transfer;
    /* dead coarse root N (kgN m-2) */
    double          deadcrootn;
    double          deadcrootn_storage;
    double          deadcrootn_transfer;
    double          cwdn;            /* coarse woody debris N (kgN m-2) */
    double          litr1n;          /* litter labile N (kgN m-2) */
    double          litr2n;          /* litter unshielded cellulose N (kgN m-2)
                                      */
    double          litr3n;          /* litter shielded cellulose N (kgN m-2) */
    double          litr4n;          /* litter lignin N (kgN m-2) */
    double          soil1n;          /* microbial recycling pool N (fast)
                                      * (kgN m-2) */
    double          soil2n;          /* microbial recycling pool N (medium)
                                      * (kgN m-2) */
    double          soil3n;          /* microbial recycling pool N (slow)
                                      * (kgN m-2) */
    double          soil4n;          /* recalcitrant SOM N (humus, slowest)
                                      * (kgN m-2) */
    double          surfn;           /* surface mineral N (kgN m-2) */
    double          sminn;           /* soil mineral N (kgN m-2) */
    double          retransn;        /* plant pool of retranslocated N (kgN m-2)
                                      */
    double          npool;           /* temporary plant N pool (kgN m-2) */
    double          nfix_src;        /* SUM of biological N fixation (kgN m-2)
                                      */
    double          ndep_src;        /* SUM of N deposition inputs (kgN m-2) */
    double          nleached_snk;    /* SUM of N leached (kgN m-2) */
    double          nvol_snk;        /* SUM of N lost to volatilization
                                      * (kgN m-2) */
    double          fire_snk;        /* SUM of N lost to fire (kgN m-2) */
} nstate_struct;

/* Daily nitrogen flux variables */
typedef struct nflux_struct
{
    /* mortality fluxes (kgN m-2 day-1) */
    double          m_leafn_to_litr1n;
    double          m_leafn_to_litr2n;
    double          m_leafn_to_litr3n;
    double          m_leafn_to_litr4n;
    double          m_frootn_to_litr1n;
    double          m_frootn_to_litr2n;
    double          m_frootn_to_litr3n;
    double          m_frootn_to_litr4n;
    double          m_leafn_storage_to_litr1n;
    double          m_frootn_storage_to_litr1n;
    double          m_livestemn_storage_to_litr1n;
    double          m_deadstemn_storage_to_litr1n;
    double          m_livecrootn_storage_to_litr1n;
    double          m_deadcrootn_storage_to_litr1n;
    double          m_leafn_transfer_to_litr1n;
    double          m_frootn_transfer_to_litr1n;
    double          m_livestemn_transfer_to_litr1n;
    double          m_deadstemn_transfer_to_litr1n;
    double          m_livecrootn_transfer_to_litr1n;
    double          m_deadcrootn_transfer_to_litr1n;
    double          m_livestemn_to_litr1n;
    double          m_livestemn_to_cwdn;
    double          m_deadstemn_to_cwdn;
    double          m_livecrootn_to_litr1n;
    double          m_livecrootn_to_cwdn;
    double          m_deadcrootn_to_cwdn;
    double          m_retransn_to_litr1n;
    /* fire fluxes (kgN m-2 day-1) */
    double          m_leafn_to_fire;
    double          m_frootn_to_fire;
    double          m_leafn_storage_to_fire;
    double          m_frootn_storage_to_fire;
    double          m_livestemn_storage_to_fire;
    double          m_deadstemn_storage_to_fire;
    double          m_livecrootn_storage_to_fire;
    double          m_deadcrootn_storage_to_fire;
    double          m_leafn_transfer_to_fire;
    double          m_frootn_transfer_to_fire;
    double          m_livestemn_transfer_to_fire;
    double          m_deadstemn_transfer_to_fire;
    double          m_livecrootn_transfer_to_fire;
    double          m_deadcrootn_transfer_to_fire;
    double          m_livestemn_to_fire;
    double          m_deadstemn_to_fire;
    double          m_livecrootn_to_fire;
    double          m_deadcrootn_to_fire;
    double          m_retransn_to_fire;
    double          m_litr1n_to_fire;
    double          m_litr2n_to_fire;
    double          m_litr3n_to_fire;
    double          m_litr4n_to_fire;
    double          m_cwdn_to_fire;
    /* phenology fluxes from transfer pools (kgN m-2 day-1) */
    double          leafn_transfer_to_leafn;
    double          frootn_transfer_to_frootn;
    double          livestemn_transfer_to_livestemn;
    double          deadstemn_transfer_to_deadstemn;
    double          livecrootn_transfer_to_livecrootn;
    double          deadcrootn_transfer_to_deadcrootn;
    /* litterfall fluxes (kgN m-2 day-1) */
    double          leafn_to_litr1n;
    double          leafn_to_litr2n;
    double          leafn_to_litr3n;
    double          leafn_to_litr4n;
    double          leafn_to_retransn;
    double          frootn_to_litr1n;
    double          frootn_to_litr2n;
    double          frootn_to_litr3n;
    double          frootn_to_litr4n;
    /* decomposition fluxes (kgN m-2 day-1) */
    double          ndep_to_sminn;
    double          nfix_to_sminn;
    /* litter and soil decomposition fluxes (kgN m-2 day-1) */
    double          cwdn_to_litr2n;
    double          cwdn_to_litr3n;
    double          cwdn_to_litr4n;
    double          litr1n_to_soil1n;
    double          sminn_to_soil1n_l1;
    double          litr2n_to_soil2n;
    double          sminn_to_soil2n_l2;
    double          litr3n_to_litr2n;
    double          litr4n_to_soil3n;
    double          sminn_to_soil3n_l4;
    double          soil1n_to_soil2n;
    double          sminn_to_soil2n_s1;
    double          soil2n_to_soil3n;
    double          sminn_to_soil3n_s2;
    double          soil3n_to_soil4n;
    double          sminn_to_soil4n_s3;
    double          soil4n_to_sminn;
    /* denitrification (volatilization) fluxes (kgN m-2 day-1) */
    double          sminn_to_nvol_l1s1;
    double          sminn_to_nvol_l2s2;
    double          sminn_to_nvol_l4s3;
    double          sminn_to_nvol_s1s2;
    double          sminn_to_nvol_s2s3;
    double          sminn_to_nvol_s3s4;
    double          sminn_to_nvol_s4;
    double          sminn_to_denitrif;
    /* leaching flux (kgN m-2 day-1) */
    double          sminn_leached;
    /* daily allocation fluxes from (kgN m-2 day-1) */
    double          retransn_to_npool;
    double          sminn_to_npool;
    double          npool_to_leafn;
    double          npool_to_leafn_storage;
    double          npool_to_frootn;
    double          npool_to_frootn_storage;
    double          npool_to_livestemn;
    double          npool_to_livestemn_storage;
    double          npool_to_deadstemn;
    double          npool_to_deadstemn_storage;
    double          npool_to_livecrootn;
    double          npool_to_livecrootn_storage;
    double          npool_to_deadcrootn;
    double          npool_to_deadcrootn_storage;
    /* annual turnover of storage to transfer pools (kgN m-2 day-1) */
    double          leafn_storage_to_leafn_transfer;
    double          frootn_storage_to_frootn_transfer;
    double          livestemn_storage_to_livestemn_transfer;
    double          deadstemn_storage_to_deadstemn_transfer;
    double          livecrootn_storage_to_livecrootn_transfer;
    double          deadcrootn_storage_to_deadcrootn_transfer;
    /* turnover of live wood to dead wood, with retranslocation (kgN m-2 day-1)
     */
    double          livestemn_to_deadstemn;
    double          livestemn_to_retransn;
    double          livecrootn_to_deadcrootn;
    double          livecrootn_to_retransn;
} nflux_struct;

/* Temporary nitrogen variables for reconciliation of decomposition
 * immobilization fluxes and plant growth N demands */
typedef struct ntemp_struct
{
    double          surfn0;             /* surface N of previous time step
                                         * (kg N m-2) */
    double          sminn0;             /* subsurface N of previous time step
                                         * (kg N m-2) */
    double          mineralized;        /* N mineralization (kgN m-2 day-1) */
    double          potential_immob;    /* potential N immobilization (kgN m-2)
                                         */
    double          plitr1c_loss;       /* potential loss from litter labile
                                         * pool (kgN m-2 day-1) */
    double          pmnf_l1s1;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          plitr2c_loss;       /* potential loss from litter unshielded
                                         * pool (kgN m-2 day-1) */
    double          pmnf_l2s2;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          plitr4c_loss;       /* potential loss from litter lignin
                                         * pool (kgN m-2 day-1) */
    double          pmnf_l4s3;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          psoil1c_loss;       /* potential loss from fast soil pool
                                         * (kgN m-2 day-1) */
    double          pmnf_s1s2;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          psoil2c_loss;       /* potential loss from medium soil pool
                                         * (kgN m-2 day-1) */
    double          pmnf_s2s3;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          psoil3c_loss;       /* potential loss from slow soil pool
                                         * (kgN m-2 day-1) */
    double          pmnf_s3s4;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    double          psoil4c_loss;       /* potential loss from slowest soil pool
                                         * (kgN m-2 day-1) */
    double          kl4;                /* decomposition rate of lignin litter
                                         * pool (day -1) */
} ntemp_struct;

/* Ecophysiological variables */
typedef struct epvar_struct
{
    double          bg_leafc_litfall_rate;  /* rate leaf litfall (kgC m-2 s-1)
                                             */
    double          bg_frootc_litfall_rate; /* rate froot litfall (kgC m-2 s-1)
                                             */
    double          livestemc_turnover_rate;/* rate livestem turnover
                                             * (kgC m-2 s-1) */
    double          livecrootc_turnover_rate; /* rate livecroot turnover
                                             * (kgC m-2 s-1) */
    double          dsr;                    /* number of days since rain */
    double          sun_proj_sla;           /* sunlit projected SLA (m2 kgC-1)
                                             */
    double          shade_proj_sla;         /* shaded projected SLA (m2 kgC-1)
                                             */
    double          psi;                    /* water potential of soil and
                                             * leaves (MPa) */
    double          dlmr_area_sun;          /* sunlit leaf MR (umolC/m2
                                             * projected leaf area s-1) */
    double          dlmr_area_shade;        /* shaded leaf MR (umolC/m2
                                             * projected leaf area s-1) */
    double          gl_t_wv_sun;            /* leaf-scale conductance to
                                             * transpired water (m s-1) */
    double          gl_t_wv_shade;          /* leaf-scale conductance to
                                             * transpired water (m s-1) */
    double          assim_sun;              /* sunlit assimilation per unit pLAI
                                             * (umol m-2 s-1) */
    double          assim_shade;            /* shaded assimilation per unit pLAI
                                             * (umol m-2 s-1) */
    double          t_scalar;               /* decomp temperature scalar (-) */
    double          w_scalar;               /* decomp water scalar (-) */
    double          rate_scalar;            /* decomp combined scalar (-) */
    double          daily_gross_nimmob;     /* daily gross N immobilization
                                             * (kgN m-2 d-1) */
    double          daily_net_nmin;         /* daily net N mineralization
                                             * (kgN m-2 d-1) */
    double          fpi;                    /* fraction of potential
                                             * immobilization (-) */
    double          m_tmin;                 /* freezing night temperature
                                             * multiplier (-) */
    double          m_psi;                  /* water potential multiplier (-) */
    double          m_co2;                  /* atmospheric (CO2) multiplier (-)
                                             */
    double          m_ppfd_sun;             /* PAR flux density multiplier (-)
                                             */
    double          m_ppfd_shade;           /* PAR flux density multiplier (-)
                                             */
    double          m_vpd;                  /* vapor pressure deficit multiplier
                                             * (-) */
    double          m_final_sun;            /* product of all other multipliers
                                             * (-) */
    double          m_final_shade;          /* product of all other multipliers
                                             * (-) */
    double          ytd_maxplai;            /* year-to-date maximum projected
                                             * LAI (-) */
    int             dormant_flag;           /* dormancy flag */
    double          days_active;            /* number of days since last
                                             * dormancy */
    int             onset_flag;             /* onset flag */
    int             onset_counter;          /* onset days counter */
    int             onset_gddflag;          /* onset flag for growing degree
                                             * day sum */
    double          onset_fdd;              /* onset freezing degree days
                                             * counter */
    double          onset_gdd;              /* onset growing degree days */
    double          onset_swi;              /* onset soil water index */
    int             offset_flag;            /* offset flag */
    int             offset_counter;         /* offset days counter */
    double          offset_fdd;             /* offset freezing degree days
                                             * counter */
    double          offset_swi;             /* offset soil water index */
    double          annavg_t2m;             /* annual average 2m air
                                             * temperature (K) */
    double          gpp;                    /* GPP flux before downregulation
                                             * (gC m-2 s-1) */
    double          prev_leafc_to_litter;   /* previous timestep leaf C
                                             * litterfall flux (gC m-2 s-1) */
    double          prev_frootc_to_litter;  /* previous timestep froot C
                                             * litterfall flux (gC m-2 s-1) */
    double          old_c_balance;          /* previous timestep C balance
                                             * (kgC m-2 day-1) */
    double          old_n_balance;          /* previous timestep N balance
                                             * (kgN m-2 day-1) */
    double          dayl;                   /* daylength (s) */
    double          prev_dayl;              /* previous day daylength (s) */
} epvar_struct;

/* Structure for the photosynthesis routine */
typedef struct psn_struct
{
    int             c3;       /* set to 1 for C3 model, 0 for C4 model */
    double          pa;       /* atmospheric pressure (Pa) */
    double          co2;      /* atmospheric CO2 concentration (ppm) */
    double          t;        /* temperature (deg C) */
    double          lnc;      /* leaf N per unit sunlit leaf area (kg Nleaf m-2)
                               */
    double          flnr;     /* fract. of leaf N in Rubisco (kg NRub/kg Nleaf)
                               */
    double          ppfd;     /* PAR flux per unit sunlit leaf area
                               * (umol m-2 s-1) */
    double          g;        /* conductance to CO2 (umol m-2 s-1 Pa-1) */
    double          dlmr;     /* day leaf maintenance respiration, projected
                               * area basis (umol m-2 s-1) */
    double          Ci;       /* intercellular CO2 concentration (Pa) */
    double          O2;       /* atmospheric O2 concentration (Pa) */
    double          Ca;       /* atmospheric CO2 concentration (Pa) */
    double          gamma;    /* CO2 compensation point, no Rd (Pa) */
    double          Kc;       /* MM constant carboxylation (Pa) */
    double          Ko;       /* MM constant oxygenation (Pa) */
    double          Vmax;     /* max rate carboxylation (umol m-2 s-1) */
    double          Jmax;     /* max rate electron transport (umol m-2 s-1) */
    double          J;        /* rate of RuBP regeneration (umol m-2 s-1) */
    double          Av;       /* carboxylation limited assimilation
                               * (umol m-2 s-1) */
    double          Aj;       /* RuBP regen limited assimilation (umol m-2 s-1)
                               */
    double          A;        /* final assimilation rate (umol m-2 s-1) */
} psn_struct;

/* CN summary structure */
typedef struct summary_struct
{
    double          daily_npp;         /* = GPP - Rmaint - Rgrowth
                                        * (kgC m-2 day-1) */
    double          daily_nep;         /* = NPP - Rheterotroph (kgC m-2 day-1)
                                        */
    double          daily_nee;         /* = NEP - fire losses (kgC m-2 day-1) */
    double          daily_gpp;         /* gross PSN source (kgC m-2 day-1) */
    double          daily_mr;          /* maintenance respiration
                                        * kgC m-2 day-1) */
    double          daily_gr;          /* growth respiration (kgC m-2 day-1) */
    double          daily_hr;          /* heterotrophic respiration
                                        * kgC m-2 day-1) */
    double          daily_fire;        /* fire losses (kgC m-2 day-1) */
    double          daily_litfallc;    /* total litterfall (kgC m-2 day-1) */
    /* summed over entire simulation (kgC m-2) */
    double          cum_npp;
    double          cum_nep;
    double          cum_nee;
    double          cum_gpp;
    double          cum_mr;
    double          cum_gr;
    double          cum_hr;
    double          cum_fire;
    double          vegc;              /* total vegetation C (kgC m-2) */
    double          litrc;             /* total litter C (kgC m-2) */
    double          soilc;             /* total soil C (kgC m-2) */
    double          totalc;            /* total of vegc, litrc, and soilc
                                        * (kgC m-2) */
} summary_struct;

/* Solute transport structure */
typedef struct solute_struct
{
    double          conc_surf;            /* surface pool concentration
                                           * (kg kgH2O-1) */
    double          conc_subsurf;         /* subsurface pool concentration
                                           * (kg kgH2O-1) */
    double          infilflux;            /* solute infiltration flux
                                           * (kg m-2 s-1) */
    double          snksrc;               /* subsurface sink/source term
                                           * (kg m-2 s-1) */
    double          ovlflux[NUM_EDGE];    /* overland solute flux (kg s-1) */
    double          subflux[NUM_EDGE];    /* subsurface solute flux (kg s-1) */
} solute_struct;
#endif

/* Element structure */
typedef struct elem_struct
{
    int             node[NUM_EDGE];    /* nodes of triangular element
                                        * (counterclockwise) */
    int             nabr[NUM_EDGE];    /* neighbor elements (neighbor i shares
                                        * edge i (0: on boundary) */
    int             ind;               /* index */
    attrib_struct   attrib;
    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
#ifdef _FBR_
    geol_struct     geol;
#endif
    epconst_struct  epc;
    ic_struct       ic;
    bc_struct       bc;
#ifdef _FBR_
    bc_struct       fbr_bc;
#endif
    wstate_struct   ws;
    wstate_struct   ws0;
    wflux_struct    wf;
    estate_struct   es;
    eflux_struct    ef;
    pstate_struct   ps;
#ifdef _DAILY_
    daily_struct    daily;
#endif
#ifdef _CYCLES_
    cropmgmt_struct cropmgmt;
    comm_struct     comm;
    residue_struct  residue;
    soilc_struct    soilc;
    weather_struct  weather;
    snow_struct     snow;
    solute_struct   NO3sol;
    solute_struct   NH4sol;
    cyclesic_struct cycles_restart;
#endif
#ifdef _BGC_
    bgcic_struct    restart_input;
    bgcic_struct    restart_output;
    cstate_struct   cs;
    cflux_struct    cf;
    nstate_struct   ns;
    nflux_struct    nf;
    psn_struct      psn_sun;
    psn_struct      psn_shade;
    ntemp_struct    nt;
    summary_struct  summary;
    epvar_struct    epv;
    solute_struct   nsol;
    spinup_struct   spinup;
#endif
} elem_struct;
#endif
