#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

/* Input file names */
typedef struct filename_struct
{
    char            riv[MAXSTRING];         /* river input file */
    char            mesh[MAXSTRING];        /* mesh structure file */
    char            att[MAXSTRING];         /* attribute file */
    char            soil[MAXSTRING];        /* soil property file */
    char            lc[MAXSTRING];          /* land cover property file */
    char            meteo[MAXSTRING];       /* meteorological forcing file */
    char            lai[MAXSTRING];         /* lai forcing file */
    char            bc[MAXSTRING];          /* boundary condition file */
    char            para[MAXSTRING];        /* control parameter file */
    char            calib[MAXSTRING];       /* calibration file */
    char            ic[MAXSTRING];          /* initial condition file */
#if defined(_BGC_)
    char            bgc[MAXSTRING];         /* bgc module control file */
    char            co2[MAXSTRING];         /* CO2 forcing file */
    char            ndep[MAXSTRING];        /* nitrogen deposition forcing file
                                             */
    char            bgcic[MAXSTRING];       /* bgc module initial condition file
                                             */
#endif
#if defined(_CYCLES_)
    char            cycles[MAXSTRING];
    char            soilinit[MAXSTRING];
    char            crop[MAXSTRING];
#endif
#if defined(_CYCLES_OBSOLETE_)
    char            op[MAXOP][MAXSTRING];
    char            cyclesic[MAXSTRING];
#endif
#if defined(_FBR_)
    char            geol[MAXSTRING];        /* geology property file */
    char            bedrock[MAXSTRING];     /* bedrock elevation file */
#endif
#if defined(_NOAH_)
    char            lsm[MAXSTRING];         /* land surface module control file
                                             */
    char            rad[MAXSTRING];         /* radiation forcing file */
    char            ice[MAXSTRING];         /* glacier ice file */
#endif
#if defined(_RT_)
    char            cdbs[MAXSTRING];        /* chemistry database file */
    char            chem[MAXSTRING];        /* chemistry control file */
    char            cini[MAXSTRING];        /* cini file */
    char            prep[MAXSTRING];        /* precipitation concentration time
                                             * series file */
    char            rtic[MAXSTRING];        /* chemistry restart file */
#endif
} filename_struct;

/* River input structure */
typedef struct rivtbl_struct
{
    int            *fromnode;               /* upstream node id */
    int            *tonode;                 /* downstream node id */
    int            *down;                   /* downstream channel id */
    int            *leftele;                /* left bank id */
    int            *rightele;               /* right bank id */
    int            *shp;                    /* river shape type */
    int            *matl;                   /* material type */
    int            *bc;                     /* boundary condition type */
    int            *rsvr;                   /* reservoir type */
} rivtbl_struct;

/* River shape parameters */
typedef struct shptbl_struct
{
    int             number;                 /* number of shape types */
    double         *depth;                  /* river channel depth */
    int            *intrpl_ord;             /* interpolation order (shape of
                                             * channel):
                                             * 1: rectangle
                                             * 2: triangle
                                             * 3: quadratic
                                             * 4: cubic */
    double         *coeff;                  /* width coefficient */
} shptbl_struct;

/* River channel material parameters */
typedef struct matltbl_struct
{
    int             number;                 /* number of bank/bed material
                                             * types */
    double         *rough;                  /* river channel roughness (s m-1/3)
                                             */
    double         *cwr;                    /* discharge coefficient (-) */
    double         *ksath;                  /* bank hydraulic conductivity
                                             * (m s-1) */
    double         *ksatv;                  /* bed hydraulic conductivity
                                             * (m s-1) */
    double         *bedthick;               /* bed thickness (m) */
} matltbl_struct;

/* Mesh structure */
typedef struct meshtbl_struct
{
    int             numnode;                /* number of nodes */
    int           **node;                   /* nodes of grids */
    int           **nabr;                   /* neighbors */
    double         *x;                      /* x of node (m) */
    double         *y;                      /* y of node (m) */
    double         *zmin;                   /* soil bottom elevation of node (m)
                                             */
    double         *zmax;                   /* surface elevation of node (m) */
#if defined(_FBR_)
    double         *zbed;                   /* impermeable bedrock elevation (m)
                                             */
#endif
} meshtbl_struct;

/* Element attribute */
typedef struct atttbl_struct
{
    int            *soil;                   /* soil type */
    int            *geol;                   /* geology type */
    int            *lc;                     /* land cover type */
    int           **bc;                     /* boundary condition type */
    int            *meteo;                  /* meteorological forcing type */
    int            *lai;                    /* leaf area index forcing type
                                             * 0: use climatological values;
                                             * else: use forcing file */
    int            *source;                 /* source forcing type */
#if defined(_FBR_)
    int           **fbr_bc;                 /* boundary condition type for
                                             * fractured bedrock layer */
#endif
#if defined(_RT_)
    int            *prcpc;                  /* precipitation concentration type
                                             */
    int           **chem_ic;                /* chemical concentration type */
#endif
} atttbl_struct;

/* Soil parameter */
typedef struct soiltbl_struct
{
    int             number;                 /* number of soil types */
    double         *silt;                   /* silt percentage (%) */
    double         *clay;                   /* clay percentage (%) */
    double         *om;                     /* organic matter percentage (%) */
    double         *bd;                     /* bulk density (g cm-3) */
    double         *kinfv;                  /* saturated infiltration
                                             * conductivity (m s-1) */
    double         *ksatv;                  /* vertical saturated hydraulic
                                             * conductivity (m s-1) */
    double         *ksath;                  /* horizontal saturated hydraulic
                                             * conductivity (m s-1) */
    double         *smcmax;                 /* maximum soil moisture content
                                             * (m3 m-3) */
    double         *smcmin;                 /* residual soil moisture content
                                             * (m3 m-3) */
    double         *smcwlt;                 /* wilting point (m3 m-3) */
    double         *smcref;                 /* soil moisture threshold where
                                             * transpiration begins to stress
                                             * (m3 m-3) */
    double         *qtz;                    /* soil quartz content (-) */
    double         *alpha;                  /* alpha from van Genuchten eqn
                                             * (m-1) */
    double         *beta;                   /* beta (n) from van Genuchten eqn
                                             * (-) */
    double         *areafh;                 /* macropore area fraction on a
                                             * horizontal cross-section (m2 m-2)
                                             */
    double         *areafv;                 /* macropore area fraction on a
                                             * vertical cross-section (m2 m-2)*/
    double         *dmac;                   /* macropore depth (m) */
    double          dinf;                   /* depth from ground surface across
                                             * which head gradient is calculated
                                             * for infiltration (m) */
    double          kmacv_ro;               /* ratio between vertical macropore
                                             * hydraulic conductivity and
                                             * vertical saturated infiltration
                                             * hydraulic conductivity (-) */
    double          kmach_ro;               /* ratio between horizontal
                                             * macropore hydraulic conductivity
                                             * and horizontal saturated
                                             * hydraulic conductivity (-) */
#if defined(_CYCLES_)
    int            *nlayers;
    double        **clay_layer;
    double        **sand_layer;
    double        **iom_layer;
    double        **bd_layer;
    double        **no3;
    double        **nh4;
    double        **b;
    double        **air_entry_pot;
    double        **fc;
    double        **pwp;
#endif
} soiltbl_struct;

/* Geology parameter */
typedef struct geoltbl_struct
{
    int             number;                 /* number of types */
    double         *ksath;                  /* horizontal saturated hydraulic
                                             * conductivity (m s-1) */
    double         *ksatv;                  /* vertical saturated hydraulic
                                             * conductivity (m s-1) */
    double         *smcmax;                 /* maximum moisture content (m3 m-3)
                                             */
    double         *smcmin;                 /* residual moisture content
                                             * (m3 m-3) */
    double         *alpha;                  /* alpha from van Genuchten eqn
                                             * (m-1) */
    double         *beta;                   /* beta (n) from van Genuchten eqn
                                             * (-) */
    double         *areafh;                 /* macropore area fraction on a
                                             * horizontal cross-section (m2 m-2)
                                             */
    double         *areafv;                 /* macropore area fraction on a
                                             * vertical cross-section (m2 m-2)*/
    double         *dmac;                   /* macropore depth (m) */
    double          kmacv_ro;               /* ratio between vertical macropore
                                             * hydraulic conductivity and
                                             * vertical saturated infiltration
                                             * hydraulic conductivity (-) */
    double          kmach_ro;               /* ratio between horizontal
                                             * macropore hydraulic conductivity
                                             * and horizontal saturated
                                             * hydraulic conductivity (-) */
} geoltbl_struct;

/* Land cover parameters */
typedef struct lctbl_struct
{
    int             number;                 /* number of land cover types */
    double         *laimax;                 /* maximum LAI across all seasons
                                             * for a vegetation type (m2 m-2) */
    double         *laimin;                 /* minimum LAI across all seasons
                                             * for a vegetation type (m2 m-2) */
    double         *vegfrac;                /* areal fractional coverage of
                                             * green vegetation (0.0-1.0) (-) */
    double         *albedomin;              /* minimum background albedo (-) */
    double         *albedomax;              /* maximum background albedo (-) */
    double         *emissmin;               /* minimum emissivity (-) */
    double         *emissmax;               /* maximum emissivity (-) */
    double         *z0min;                  /* minimum roughness length (m) */
    double         *z0max;                  /* maximum roughness length (m) */
    double         *hs;                     /* parameter used in vapor pressure
                                             * deficit function (-) */
    double         *snup;                   /* threshold snow depth (in water
                                             * equivalent) that implies 100%
                                             * snow cover (m) */
    double         *rgl;                    /* reference incoming solar flux for
                                             * photosynthetically active canopy
                                             * (W m-2) */
    double         *rsmin;                  /* minimum canopy resistance (s m-1)
                                             */
    double         *rough;                  /* surface roughness (Manning's n)
                                             * (s m-1/3) */
    double         *rzd;                    /* rooting depth (m) */
    double          rsmax;                  /* cuticular resistance (s m-1) */
    double          cfactr;                 /* parameter used in the canopy
                                             * interception calculation (-) */
    double          topt;                   /* optimum transpiration air
                                             * temperature (K) */
} lctbl_struct;

/* Time series data structure */
typedef struct tsdata_struct
{
    int             length;                 /* length of time series */
    int            *ftime;                  /* forcing time */
    double        **data;                   /* forcing values at forcing time */
    double         *value;                  /* forcing values at model time t */
    union
    {
        int             bc_type;            /* boundary condition type:
                                             * 1 = Dirichlet, 2 = Neumann */
        double          zlvl_wind;          /* height above ground of wind
                                             * observations (m) */
    };
} tsdata_struct;

/* Forcing structure */
typedef struct forc_struct
{
    int             nbc;                    /* number of boundary condition
                                             * series */
    tsdata_struct  *bc;                     /* boundary condition time series */
    int             nmeteo;                 /* number of meteorological forcing
                                             * series */
    tsdata_struct  *meteo;                  /* meteorological forcing series */
    int             nlai;                   /* number of lai series */
    tsdata_struct  *lai;                    /* lai forcing series */
    int             nsource;                /* number of source forcing series*/
    tsdata_struct  *source;                 /* source forcing series */
    int             nriverbc;               /* number of river boundary
                                             * conditions */
    tsdata_struct  *riverbc;                /* river boundary condition series*/
#if defined(_BGC_)
    int             nco2;
    tsdata_struct  *co2;                    /* CO2 forcing series */
    int             nndep;
    tsdata_struct  *ndep;                   /* nitrogen deposition forcing
                                             * series */
#endif
#if defined(_NOAH_)
    int             nrad;                   /* number of radiation forcing
                                             * series */
    tsdata_struct  *rad;                    /* radiation forcing series */
#endif
#if defined(_RT_)
    int             PrpFlg;                 /* flag that indicates how
                                             * precipitation is specified */
    int             nprcpc;                 /* number of precipitation
                                             * concentration time series */
    tsdata_struct  *prcpc;                  /* concentration in precipitation */
#endif
} forc_struct;

#if defined(_NOAH_)
/* Land surface parameters */
typedef struct noahtbl_struct
{
    double          sbeta;                  /* parameter used to calculate
                                             * vegetation effect on soil heat
                                             * (-) */
    double          fxexp;                  /* soil evaporation exponent used in
                                             * direct evaporation (-) */
    double          csoil;                  /* soil heat capacity (J m-3 K-1) */
    double          salp;                   /* shape parameter of distribution
                                             * function of snow cover (-) */
    double          frzk;                   /* frozen ground parameter (-) */
    double          zbot;                   /* depth of lower boundary soil
                                             * temperature (m) */
    double          tbot;                   /* bottom soil temperature (local
                                             * yearly-mean surface air
                                             * temperature) (K) */
    double          czil;                   /* Zilitinkevich constant (-) */
    double          lvcoef;                 /* parameter controls surface snow
                                             * albedo in the presence of
                                             * snowcover (-) */
} noahtbl_struct;
#endif

#if defined(_BGC_)
/* Ecophysiological parameters */
typedef struct epctbl_struct
{
    int            *woody;                  /* flag: 1 = woody, 0 = non-woody */
    int            *evergreen;              /* flag: 1 = evergreen,
                                             * 0 = deciduous */
    int            *c3_flag;                /* flag: 1 = C3,  0 = C4 */
    int            *phenology_flag;         /* flag: 1 = phenology model,
                                             * 0 = user defined */
    int            *onday;                  /* day of year when leaves on */
    int            *offday;                 /* day of year when leaves off */
    int            *transfer_days;          /* growth period for transfer (day)
                                             */
    int            *litfall_days;           /* growth period for litfall (day)
                                             */
    double         *leaf_turnover;          /* annual leaf turnover fraction
                                             * (yr-1) */
    double         *froot_turnover;         /* annual fine root turnover
                                             * fraction (yr-1) */
    double         *livewood_turnover;      /* annual live wood turnover
                                             * fraction (yr-1) */
    double         *daily_mortality_turnover;/* daily mortality turnover
                                             * (day-1) */
    double         *daily_fire_turnover;    /* daily fire turnover (day-1) */
    double         *alloc_frootc_leafc;     /* new fine root C to new leaf C (-)
                                             */
    double         *alloc_newstemc_newleafc;/* new stem C to new leaf C (-) */
    double         *alloc_newlivewoodc_newwoodc;/* new livewood C:new wood C
                                             * (-) */
    double         *alloc_crootc_stemc;     /* new live croot C to new live stem
                                             * C (-) */
    double         *alloc_prop_curgrowth;   /* daily allocation to current
                                             * growth (-) */
    double         *avg_proj_sla;           /* canopy average projected SLA
                                             * (m2 kgC-1) */
    double         *sla_ratio;              /* ratio of shaded to sunlit
                                             * projected SLA (-) */
    double         *lai_ratio;              /* ratio of (all-sided LA /
                                             * one-sided LA) (-) */
    double         *ext_coef;               /* canopy light extinction
                                             * coefficient (-) */
    double         *flnr;                   /* leaf N in Rubisco
                                             * (kgNRub kgNleaf-1) */
    double         *psi_open;               /* psi at start of conductance
                                             * reduction (MPa) */
    double         *psi_close;              /* psi at complete conductance
                                             * reduction (MPa) */
    double         *vpd_open;               /* vpd at start of conductance
                                             * reduction (Pa) */
    double         *vpd_close;              /* vpd at complete conductance
                                             * reduction (Pa) */
    double         *froot_cn;               /* C:N for fine roots (kgC kgN-1) */
    double         *leaf_cn;                /* C:N for leaves (kgC kgN-1) */
    double         *livewood_cn;            /* C:N for live wood (kgC kgN-1) */
    double         *deadwood_cn;            /* C:N for dead wood (kgC kgN-1) */
    double         *leaflitr_cn;            /* constant C:N for leaf litter
                                             * (kgC kgN-1) */
    double         *leaflitr_flab;          /* leaf litter labile fraction (-)*/
    double         *leaflitr_fucel;         /* leaf litter unshielded cellulose
                                             * fraction (-) */
    double         *leaflitr_fscel;         /* leaf litter shielded cellulose
                                             * fraction (-) */
    double         *leaflitr_flig;          /* leaf litter lignin fraction (-)*/
    double         *frootlitr_flab;         /* fine root litter labile fraction
                                             * (-) */
    double         *frootlitr_fucel;        /* fine root litter unshielded
                                             * cellulose fraction (-) */
    double         *frootlitr_fscel;        /* fine root litter shielded
                                             * cellulose fraction (-) */
    double         *frootlitr_flig;         /* fine root litter lignin fraction
                                             * (-) */
    double         *deadwood_fucel;         /* dead wood unshielded cellulose
                                             * fraction (-) */
    double         *deadwood_fscel;         /* dead wood shielded cellulose
                                             * fraction (-) */
    double         *deadwood_flig;          /* dead wood lignin fraction (-) */
} epctbl_struct;
#endif

#if defined(_CYCLES_)
typedef struct agtbl_struct
{
    int            *oper;
    int             noper;
    char            oper_filen[MAXOP][MAXSTRING];
} agtbl_struct;

typedef struct crop_epc_struct
{
    char            name[MAXSTRING];        /* name of crop (-) */
    double          flower_thermal_time;    /* thermal time to flowering
                                             * (degree C day) */
    double          mat_thermal_time;       /* thermal time to maturity
                                             * (degree C day) */
    double          max_soil_cover;         /* maximum crop cover (-) */
    double          max_rooting_depth;      /* maximum rooting depth (m) */
    double          frac_residue_stand;     /* fraction of aboveground residues
                                             * in standing position (-) */
    double          frac_residue_removed;   /* fraction of non-grain crop
                                             * biomass removed with harvested
                                             * grain; or fraction of harvestable
                                             * aboveground biomass removed with
                                             * harvested forages or grazed (-)*/
    double          clip_biomass_thld_upper;/* aboveground plant biomass
                                             * threshold that triggers clipping
                                             * event or forage harvest (Mg ha-1)
                                             */
    double          clip_biomass_thld_lower;/* aboveground plant biomass
                                             * threshold that remains
                                             * un-harvestable during clipping
                                             * events (Mg ha-1) */
    double          clip_timing;            /* fraction of thermal time to crop
                                             * maturity that triggers clipping
                                             * event or grain harvest (-) */
    int             kill_after_harvest;     /* kill after harvest flag */
    int             clip_density;           /* destiny of biomass cut by
                                             * clipping events (-) */
    double          tranp_tmp_min;          /* air temperature below which
                                             * transpiration ceases (degree C)*/
    double          tranp_tmp_thld;         /* threshold air temperature for
                                             * transpiration calculation
                                             * (degree C) */
    double          cold_damage_tmp_min;    /* minimum temperature in cold
                                             * damage factor (degree C) */
    double          cold_damage_tmp_thld;   /* threshold temperature in cold
                                             * damage factor (degree C) */
    double          tmp_base;               /* base temperature for phenological
                                             * development (degree C) */
    double          tmp_opt;                /* optimum temperature for
                                             * phenological development
                                             * (degree C) */
    double          tmp_max;                /* maximum temperature for
                                             * phenological development
                                             * (degree C) */
    double          shoot_partn_init;       /* fraction of growth partitioned to
                                             * shoot biomass at emergence (-) */
    double          shoot_partn_final;      /* fraction of growth partitioned to
                                             * shoot biomass at maturity (-) */
    double          rad_use_eff;            /* radiation use efficiency (g MJ-1)
                                             */
    double          transp_use_eff;         /* transpiration use efficiency at
                                             * 1 kPa VPD (g kg-1) */
    double          hi_max;                 /* maximum harvest index (-) */
    double          hi_opt;                 /* optimum harvest index (-) */
    double          hi_min;                 /* minimum harvest index (-) */
    double          emergen_thermal_time;   /* thermal time to emergence
                                             * (degree C) */
    double          n_conc_max;             /* maximum N concentration (g g-1)*/
    double          n_dil_slope;            /* N dilution curve slope parameter
                                             * (-) */
    double          kc;                     /* transpiration coefficient (-) */
    int             annual;                 /* annual/perennial flag (-) */
    int             legume;                 /* legume flag (-) */
    int             c3;                     /* C3/C4 flag (-) */
    double          lwp_stress_onset;       /* leaf water potential for onset of
                                             * stress (J kg-1) */
    double          lwp_wilting_point;      /* leaf water potential at wilting
                                             * point (J kg-1) */
    double          tranp_max;              /* maximum transpiration rate
                                             * (mm day-1) */
} crop_epc_struct;

typedef struct realized_crop_struct
{
    char            name[MAXSTRING];        /* name of crop */
    int             plant_ymd;              /* year and date of planting */
    double          biomass;                /* total biomass at harvest
                                             * (Mg ha-1) */
    double          root;                   /* root biomass at harvest (Mg ha-1)
                                             */
    double          grain_yield;            /* grain yield at harvest (Mg ha-1)
                                             */
    double          forage_yield;           /* forage yield at harvest (Mg ha-1)
                                             */
    double          residue;                /* aboveground residue biomass left
                                             * in field at harvest (Mg ha-1) */
    double          harvest_index;          /* fraction of aboveground biomass
                                             * harvested as grain (-) */
    double          n_total;                /* total biomass nitrogen at harvest
                                             * (Mg ha-1) */
    double          n_root;                 /* root biomass nitrogen at harvest
                                             * (Mg ha-1) */
    double          n_yield_grain;          /* grain nitrogen content at harvest
                                             * (Mg ha-1) */
    double          n_yield_forage;         /* forage or removed residue N
                                             * content at harvest (Mg ha-1) */
    double          n_stress_cum;           /* cumulative N stress over duration
                                             * of crop growth (-) */
    double          n_in_harvest;           /* N content in removed biomass
                                             * (kg ha-1) */
    double          n_in_residue;           /* N content in residue left in
                                             * field (kg ha-1) */
    double          n_conc_forage;          /* N concentration in forage or
                                             * removed residues (%) */
    double          n_auto_added;           /* cumulative N added in auto
                                             * fertilization (Mg ha-1) */
    double          n_fix;                  /* cumulative N fixation by legume
                                             * (Mg ha-1) */
    double          transp;                 /* crop transpiration (mm) */
    double          transp_pot;             /* potential crop transpiration (mm)
                                             */
    double          soil_evap;              /* soil evaporation (mm) */
} realized_crop_struct;

typedef struct crop_struct
{
    crop_epc_struct epc;
    int             stage_growth;           /* phenological stage */
    int             auto_irrig;             /* auto irrigation index */
    int             auto_fert;              /* auto fertilization flag */
    int             plant_ymd;              /* year and date of planting */
    double          plant_density;          /* see planting struct */
    int             clip_start;             /* see planting struct */
    int             clip_end;               /* see planting struct */
    /* State Variables */
    double          thermal_time_daily;     /* daily thermal time (degree C) */
    double          thermal_time_cum;       /* cumulative thermal time
                                             * (degree C) */
    double          rad_intcp;              /* solar radiation intercepted by
                                             * green leaves (-) */
    double          rad_intcp_brown;        /* solar radiation intercepted
                                             * by brown leaves (-) */
    double          rad_intcp_nc;           /* solar radiation intercepted by
                                             * green leaves without considering
                                             * competition (-) */
    double          biomass;                /* total biomass (Mg ha-1) */
    double          shoot;                  /* shoot biomass (Mg ha-1) */
    double          root;                   /* root biomass (Mg ha-1) */
    double          rhizo;                  /* rhizo biomass (Mg ha-1) */
    double          shoot_growth;           /* daily shoot growth (Mg ha-1) */
    double          root_growth;            /* daily root growth (Mg ha-1) */
    double          rhizo_deposit;          /* daily rhizo deposition (Mg ha-1)
                                             */
    double          shoot_growth_unstr;     /* unstressed daily shoot
                                             * growth (Mg ha-1) */
    double          root_growth_unstr;      /* unstressed daily root growth
                                             * (Mg ha-1) */
    double          shoot_post_flower;      /* shoot biomas cumulated after
                                             * flowering (Mg ha-1) */
    double          rooting_depth;          /* rooting depth (m) */
    double          transp;                 /* daily transpiration (mm day-1)*/
    double          transp_pot;             /* daily potential transpiration
                                             * (mm day-1) */
    double          n_shoot;                /* shoot biomass N content (Mg ha-1)
                                             */
    double          n_root;                 /* root biomass N content (Mg ha-1)
                                             */
    double          n_rhizo;                /* rhizo N content (Mg ha-1) */
    double          n_rhizo_deposit;        /* rhizo daily N deposition
                                             * (Mg ha-1) */
    double          n_auto_added;           /* cumulative N added in auto
                                             * fertilization (Mg ha-1) */
    double          n_fix;                  /* cumulative N fixation by legume
                                             * (Mg ha-1) */
    double          water_stress;           /* daily water stress factor (-) */
    double          n_stress;               /* daily nitrogen stress factor (-)
                                             */
    double          shoot_unstr;            /* cumulative unstressed shoot
                                             * growth (Mg ha-1) */
    double          n_stress_cum;           /* cumulative N stress factor (-) */
    int             harvest_date_final;     /* final harvest day of year (-) */
    int             harvest_count;          /* total count of harvests (-) */
    realized_crop_struct rc;
} crop_struct;

#endif

#if defined(_CYCLES_OBSOLETE_)
typedef struct epconst_struct
{
    char            cropn[MAXSTRING];
    double          flowering_tt;           /* (degree C day) */
    double          maturity_tt;            /* (degree C day) */
    double          max_soil_covr;          /* (100%) */
    double          max_root_dpth;          /* (m) */
    double          frac_res_standing;      /* (100%) */
    double          frac_res_removed;       /* (100%) */
    double          clip_biomass_ut;        /* (kg m-2) */
    double          clip_biomass_lt;        /* (kg m-2) */
    double          clip_timing;            /* (100% thermal time) */
    int             clip_destiny;           /* (-) */
    double          transp_temp_min;        /* (degree C) */
    double          transp_temp_thr;        /* (degree C) */
    double          cold_damage_temp_min;   /* (degree C) */
    double          cold_damagetemp_thr;    /* (degree C) */
    double          temp_base;              /* (degree C) */
    double          temp_opt;               /* (degree C) */
    double          temp_max;               /* (degree C) */
    double          shoot_par_init;         /* (100%) */
    double          shoot_par_final;        /* (100%) */
    double          rue;                    /* kg J-1 */
    double          tue;                    /* kg kgH2O-1 at VPD = 1kPa */
    double          hi_max;                 /* (-) */
    double          hi_min;                 /* intercept harvest index (-) */
    double          hi;                     /* (-) */
    double          emergence_tt;           /* (degree C day) */
    double          n_conc_max;             /* (g g-1) */
    double          n_diln_slope;           /* (-) */
    double          kc;                     /* (-) */
    int             annual;                 /* (-) */
    int             legume;                 /* (-) */
    int             c3c4;                   /* (-) */
    double          lwp_stress_onset;       /* (J kg-1, or m2 s-2) */
    double          lwp_wlt;                /* (J kg-1, or m2 s-2) */
    double          transp_max;             /* (mm day-1) */
} epconst_struct;

typedef struct plant_struct
{
    /* Planting */
    int             year;
    int             doy;
    int             crop_id;
    int             auto_irr;
    int             auto_fert;
    double          plant_density;          /* (100%) */
    int             clip_start;             /* (day of year) */
    int             clip_end;               /* (day of year) */
    int             ai_start;
    int             ai_stop;
    double          ai_h2o_depl;            /* (100% plant avialable water
                                             * content) */
    int             ai_last_lyr;
} plant_struct;

typedef struct tillage_struct
{
    /* Tillage */
    int             year;
    int             doy;
    char            tooln[MAXSTRING];
    double          depth;                  /* (m) */
    double          sdr;                    /* (-) */
    double          mix_eff;                /* (100%) */
    double          drop_id;
#if NOT_YET_IMPLEMENTED
    double          fractionThermalTime;
    double          killEfficiency;
#endif
    int             grain_harv;
    double          forage_harv;
} tillage_struct;

typedef struct fixirr_struct
{
    /* Fixed Irrigation */
    int             year;
    int             doy;
    double          volume;                 /* (mm) */
} fixirr_struct;

typedef struct fixfert_struct
{
    /* Fixed Fertilization */
    int             year;
    int             doy;
    char            source[MAXSTRING];
    double          mass;                   /* (kg m-2) */
#if NOT_YET_IMPLEMENTED
    char            opForm[MAXSTRING];
    char            opMethod[MAXSTRING];
#endif
    int             layer;                  /* Starting from 1 */
    double          c_org;                  /* (100%) */
    double          c_cac;                  /* (100%) */
    double          n_org;                  /* (100%) */
    double          n_cac;                  /* (100%) */
    double          nh4;                    /* (100%) */
    double          no3;                    /* (100%) */
    double          p_org;                  /* (100%) */
    double          p_cac;                  /* (100%) */
    double          p_inorg;                /* (100%) */
    double          k;                      /* (100%) */
    double          s;                      /* (100%) */
} fixfert_struct;

typedef struct autoirr_struct
{
    int             crop_id;
    int             start;
    int             stop;
    double          h2o_depl;               /* (100% plant avialable water
                                             * content) */
    int             last_lyr;
} autoirr_struct;

typedef struct opertbl_struct
{
    fixfert_struct *fix_fert;
    int             nfert;
    fixirr_struct  *fix_irr;
    int             nirr;
    tillage_struct *tillage;
    int             ntill;
    plant_struct   *plant;
    int             nplant;
    autoirr_struct *auto_irr;
    int             nauto_irr;
} opertbl_struct;
#endif

#if defined(_RT_)
typedef struct chemtbl_struct
{
    char            name[MAXSTRING];        /* molecular formula or name */
    double          molar_mass;             /* (g mol-1) */
    double          molar_vol;              /* (cm3 mol-1) */
    double          charge;                 /* charge */
    double          size_fac;               /* size factor for DH equation */
    int             itype;                  /* type of primary species
                                             * 1 = primary aqueous
                                             * 2 = primary adsorption
                                             * 3 = primary cation exchange
                                             * 4 = primary mineral */
    int             mtype;                  /* type of the mass action species
                                             * 0 = immobile mass action
                                             * 1 = mobile mass action
                                             * 2 = mixed mobility mass action */
} chemtbl_struct;

typedef struct kintbl_struct
{
    int             position;               /* position of target mineral in the
                                             * array of primary species */
    char            label[MAXSTRING];       /* label of kinetic reaction */
    int             type;                   /* type of the kinetic reaction
                                             * 1: tst  2: precipitation-only
                                             * 3: dissolution-only 4: monod */
    double          rate;                   /* rate of kinetic reaction */
    double          actv;                   /* activation energy, used to
                                             * calculate kinetic rate of
                                             * reaction under different
                                             * temperatures */
    double          Keq;                    /* equilibrium constant */
    int             ndep;                   /* number of dependency */
    int             dep_index[MAXDEP];      /* position of species that kinetic
                                             * reaction depends on */
    double          dep_power[MAXDEP];      /* power of dependency */
    int             biomass_index;
    int             nmonod;
    int             monod_index[MAXDEP];
    double          monod_para[MAXDEP];     /* parameter for monod dependency */
    int             ninhib;
    int             inhib_index[MAXDEP];
    double          inhib_para[MAXDEP];     /* parameters that controls this
                                             * inhibition */
} kintbl_struct;

typedef struct rttbl_struct
{
    int             actv_mode;              /* activity coefficient mode:
                                             * 0 = unity coefficient,
                                             * 1 = DH equation */
    int             tmp_coup;               /* flag to couple soil temperature*/
    int             rel_min;                /* relative mineral flag:
                                             * 1 = total solid volume,
                                             * 0 = total pore volume */
    int             transpt_flag;           /* transport only flag:
                                             * 0 = simulate kinetic reaction,
                                             * 1 = transport only */
    double          cond;                   /* ratio between infiltration
                                             * concentration and rain water */
    int             num_stc;                /* number of total species */
    int             num_spc;                /* number of primary species */
    int             num_ssc;                /* number of secondary speices */
    int             num_sdc;                /* number of independent species */
    int             num_min;                /* number of minerals */
    int             num_ads;                /* number of adsorption species */
    int             num_cex;                /* number of cation exchange */
    int             num_mkr;                /* number of mineral kinetic
                                             * reactions */
    int             num_akr;                /* number of aqueous kinetic
                                             * reactions */
    double          diff_coef;              /* diffusion coefficient (m2 s-1) */
    double          disp_coef;              /* dispersion coefficient (m) */
    double          cementation;            /* cementation factor that
                                             * represents connectivity of pores
                                             */
    double          tmp;                    /* temperature of the moment */
    double          prcp_conc[MAXSPS];
    double          dep_mtx[MAXSPS][MAXSPS];/* dependency of secondary
                                             * species on primary species */
    double          dep_kin[MAXSPS][MAXSPS];/* dependency of kinetic species
                                             * on primary species  */
    double          conc_contrib[MAXSPS][MAXSPS];/* contribution of each species
                                             * to total concentration */
#if NOT_YET_IMPLEMENTED
    double          Totalconck[MAXSPS][MAXSPS];/* contribution of each species
                                             * to total concentration with
                                             * kinetic reaction included */
#endif
    double          keq[MAXSPS];            /* Keq's of secondary species */
    double          keq_kin[MAXSPS];      /* Keq's of kinetic species */
    double          adh;                    /* Debye Huckel parameter */
    double          bdh;                    /* Debye Huckel parameter */
    double          bdt;                    /* Debye Huckel parameter */
} rttbl_struct;

typedef struct chmictbl_struct
{
    int             nic;
    double        **conc;                   /* chemical concentration */
    double        **ssa;                    /* specific surface area */
} chmictbl_struct;
#endif

#endif
