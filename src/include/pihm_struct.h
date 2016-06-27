#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

/* 
 * Topographic parameters
 */
typedef struct topo_struct
{
    double          area;       /* area of element */
    double          x;          /* x of centroid */
    double          y;          /* y of centroid */
    double          zmin;       /* z_min of centroid */
    double          zmax;       /* z_max of centroid */
    double          zbed;
    double          edge[3];    /* edge i is from node i to node i+1 */
    double          surfx[3];
    double          surfy[3];
    double          node_zmax;
#ifdef _NOAH_
    double          slope;
    double          aspect;     /* ys: surface aspect of grid */
    double          svf;        /* ys: sky view factor */
    double          h_phi[36];  /* unobstrcted angle */
#endif
#ifdef _RT_
    double          areasub[3];
#endif
} topo_struct;

/*
 * Soil parameters
 */
typedef struct soil_struct
{
    int             type;

    double          depth;
    double          ksath;      /* horizontal geologic saturated
                                 * hydraulic conductivity */
    double          ksatv;      /* vertical geologic saturated
                                 * hydraulic conductivity */
    double          kinfv;      /* vertical surface saturated hydraulic
                                 * conductivity */
    double          dinf;       /* depth from ground surface accross
                                 * which head is calculated during
                                 *  infiltration */
    double          alpha;      /* alpha from van Genuchten eqn */
    double          beta;
    double          porosity;
    double          smcmax;
    double          smcmin;
    double          smcwlt;     /* wilting point (volumetric) */
    double          smcref;     /* soil moisture threshold where transpiration begins to stress (volumetric) */
    double          dmac;       /* macropore Depth */
    double          kmach;      /* macropore horizontal saturated
                                 * hydraulic conductivity */
    double          kmacv;      /* macropore vertical saturated
                                 * hydraulic conductivity */
    double          areafv;     /* macropore area fraction on a
                                 * vertical cross-section */
    double          areafh;     /* macropore area fraction on a
                                 * horizontal cross-section */
#ifdef _NOAH_
    double          csoil;      /* soil heat capacity (j m-3 k-1) */
    double          quartz;     /* soil quartz content */
    double          smcdry;     /* dry soil moisture threshold where direct evap frm top layer ends (volumetric) */
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
//    double          FC_WaterPotential[MAXLYR];
//    double          airEntryPotential[MAXLYR];
//    double          B_Value[MAXLYR];
//    double          M_Value[MAXLYR];
//    double          ksat[MAXLYR];

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
    double          NO3Leaching;
    double          NH4Leaching;

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
#endif
} soil_struct;

/*
 * Land cover ecophysiological paraemters
 */
typedef struct lc_struct
{
    int             type;
    double          shdfac;     /* areal fractional coverage of green vegetation (fraction= 0.0-1.0) */
    double          shdmin;     /* minimum areal fractional coverage of green vegetation (fraction= 0.0-1.0) <= shdfac */
    double          shdmax;     /* maximum areal fractional coverage of green vegetation (fraction= 0.0-1.0) >= shdfac */
    double          rzd;
    double          rsmin;      /* minimum canopy resistance */
    double          rgl;        /* reference incoming solar flux for
                                 * photosynthetically active canopy */
    double          laimin;
    double          laimax;     /* maxm. LAI accross all seasons for a
                                 * vegetation type */
    double          snup;       /* threshold snow depth (in water equivalent m) that implies 100 percent snow cover */
    double          hs;         /* parameter used in vapor pressure deficit function */
    double          topt;       /* optimum transpiration air temperature */
    double          cfactr;     /* parameter used in the canopy inteception calculation */
    double          rsmax;      /* max. stomatal resistance */
    double          emissmax;   /* maximum emissivity */
    double          emissmin;   /* minimum emmisivity */
    double          albedomax;  /* maximum background albedo */
    double          albedomin;  /* manimum background albedo */
    double          z0max;      /* maximum roughness length (m) */
    double          z0min;      /* minimum roughness length (m) */
    double          rough;      /* surface roughness of an element */
    double          cmcfactr;
    int             bare;
#ifdef _NOAH_
    int             isurban;
    int             nroot;      /* number of root layers, a function of veg type, determined in subroutine redprm. */
    double          ptu;        /* photo thermal unit (plant phenology for annuals/crops) (not yet used, but passed to redprm for future use in veg parms) */
    double          rtdis[MAXLYR];
#endif
} lc_struct;

/*
 * Land cover physical states
 */
typedef struct pstate_struct
{
    int             macpore_status;
    double          rc;         /* lcp: canopy resistance (s m-1) */
    double          pc;         /* plant coefficient (unitless fraction, 0-1) where pc*etp = actual transp */
    double          xlai;       /* lcp: leaf area index (dimensionless) */
    double          rcs;        /* incoming solar rc factor (dimensionless) */
    double          rct;        /* air temperature rc factor (dimensionless) */
    double          rcq;        /* atmos vapor pressure deficit rc factor (dimensionless) */
    double          rcsoil;     /* soil moisture rc factor (dimensionless) */
    double          alb;        /* backround snow-free surface albedo (fraction), for julian day of year (usually from temporal interpolation of monthly */
    double          albedo;
    double          zlvl;       /* lcp: height (m) above ground of atmospheric forcing
                                 * variables */
    double          zlvl_wind;  /* lcp: height (m) above ground of wind observations */
    double          sfcspd;
    double          rh;
    double          sfcprs;
#ifdef _NOAH_
    double          snoalb;     /* upper bound on maximum albedo over deep snow (e.g. from robinson and kukla, 1985, j. clim. & appl. meteor.) */
    int             nsoil;      /* number of soil layers (at least 2, and not
                                 * greater than parameter nsold set below) */

    double          sldpth[MAXLYR];     /* the thickness of each soil layer (m) */
    double          frzk;       /* frozen ground parameter */
    double          frzx;
    double          czil;       /* calculate roughness length of heat */
    double          emissi;     /* lcp: surface emissivity (between 0 and 1) */
    double          ch;         /* lcp: surface exchange coefficient for heat and moisture (m s-1); */
    double          rch;
    /*
     * note: ch is technically a conductance since it has been multiplied by wind speed. 
     */
    double          cm;         /* lcp: surface exchange coefficient for momentum (m s-1); */
    /*
     * note: cm is technically a conductance since it has been multiplied by wind speed. 
     */
    double          z0;         /* lcp: time varying roughness length (m) as function of snow depth */
    double          fcr;        /* YS: reduction of infiltration caused
                                 * by frozen ground */
    int             nmacd;
    double          salp;       /* shape parameter of distribution function of snow cover */
    double          fxexp;      /* soil evaporation exponent used in devap */
    double          sbeta;      /* parameter used to calculate vegetation effect on soil heat */
    double          lvcoef;
    double          snotime1;   /* age of the snow on the ground */
    double          ribb;       /* bulk richardson number used to limit the dew/frost */
    double          beta;       /* ratio of actual/potential evap (dimensionless) */
    double          sncovr;     /* fractional snow cover (unitless fraction, 0-1) */
    double          q1;         /* effective mixing ratio at surface (kg kg-1), used for diagnosing the mixing ratio at 2 meter for coupled model */
    double          q2;         /* mixing ratio at height zlvl above ground (kg kg-1) */
    double          cosz;       /* solar zenith angle (not used for now) */
    double          ffrozp;     /* fraction of frozen precipitation */
    double          z0brd;      /* background fixed roughness length (m) */
    double          embrd;      /* background surface emissivity (between 0 and 1) */
    double          q2sat;
    double          dqsdt2;
    int             nwtbl;
    double          sndens;
    double          snowh;      /* actual snow depth (m) */
    double          sncond;
    double          rr;         /* (-) */
    double          eta_kinematic;      /* atctual latent heat flux in kg m-2 s-1 */
    double          zbot;       /* depth (m) of lower boundary soil temperature */
    double          tbot;       /* bottom soil temperature (local yearly-mean sfc air temperature) */
    double          gwet;
    double          satdpth[MAXLYR];
#endif
} pstate_struct;

typedef struct elemforc_struct
{
    int             bc_type[3];
    int             lai_type;

    /* Forcing values */
    double         *bc[3];
    double         *meteo[NUM_METEO_VAR];
    double         *lai;
    double         *z0;
    double         *source;
    double         *meltf;
    double         *riverbc;
#ifdef _NOAH_
    double         *rad[2];
#endif
} elemforc_struct;

#ifdef _DAILY_
typedef struct daily_struct
{
    int             counter;
    int             daylight_counter;

    double          sfctmp;
    double          tday;
    double          tnight;
    double          tmax;
    double          tmin;

    double          solar;
    double          solar_total;

    double          sfcspd;

    double          sfcprs;

    double          surf;
    double          unsat;
    double          gw;

#ifdef _CYCLES_
    double          et[MAXLYR];
    double          sncovr;
#endif
    double          fluxsurf[4];
    double          fluxsub[4];
    double          infil;
    double          rechg;

    double          dayl;
    double          prev_dayl;

    double          stc[MAXLYR];
    double          sh2o[MAXLYR];
    double          smflxv[MAXLYR];
    double          smflxh[4][MAXLYR];
    double          q2d;
    double          albedo;
    double          ch;
} daily_struct;
#endif

typedef struct wstate_struct
{
    double          stage;
    double          surf;
    double          gw;
    double          unsat;
    double          sneqv;      /* liquid water-equivalent snow depth (m). note: snow density = sneqv/snowh */
    double          cmcmax;     /* maximum canopy water capacity */
    double          cmc;        /* Interception storage */
#ifdef _NOAH_
    double          smc[MAXLYR];        /* total soil moisture content (volumetric fraction) */
    double          sh2o[MAXLYR];       /* unfrozen soil moisture content (volumetric fraction). note: frozen soil moisture = smc - sh2o */
    double          soilw;      /* available soil moisture in root zone (unitless fraction between smcwlt and smcmax) */
    double          soilm;      /* total soil column moisture content (frozen+unfrozen) (m) */
#endif
} wstate_struct;

typedef struct wflux_struct
{
    double          fluxsurf[3];        /* Overland Flux */
    double          fluxsub[3]; /* Subsurface Flux */
    double          prcp;       /* Precep. on each element */
    double          netprcp;    /* Net precep. on each elment */
    double          infil;      /* Variable infiltration rate */
    double          rechg;      /* Recharge rate to GW */
    double          drip;       /* through-fall of precip and/or dew in excess of canopy water-holding capacity (m/s) */
    double          edir;
    double          ett;
    double          ec;
    double          etp;
    double          eta;
    double          edir_sfc;
    double          edir_unsat;
    double          edir_gw;
    double          ett_unsat;
    double          ett_gw;
    double          fluxriv[11];
#ifdef _NOAH_
    double          et[MAXLYR];
    double          runoff1;    /* surface runoff (m s-1), not infiltrating the surface */
    double          runoff2;    /* subsurface runoff (m s-1), drainage out bottom of last soil layer (baseflow) */
    double          runoff2_lyr[MAXLYR];
    double          runoff3;    /* numerical trunctation in excess of porosity (smcmax) for a given soil layer at the end of a time step (m s-1). note: the above runoff2 is actually the sum of runoff2 and runoff3 */
    double          smflxv[MAXLYR];
    double          smflxh[4][MAXLYR];
    double          pcpdrp;     /* combined prcp1 and drip (from cmc) that goes into the soil (m s-1) */
    double          prcprain;   /* liquid-precipitation rate (kg m-2 s-1) (not used) */
    double          dew;        /* dewfall (or frostfall for t<273.15) (m) */
    double          snomlt;     /* snow melt (m/s) (water equivalent) */
    double          esnow;      /* sublimation from (or deposition to if <0) snowpack (ms-1) */
    double          etns;
#endif
#ifdef _CYCLES_
    double          eres;
#endif
} wflux_struct;

typedef struct estate_struct
{
    double          stc[MAXLYR];        /* soil temp (k) */
    double          t1;         /* ground/canopy/snowpack) effective skin temperature (k) */
    double          th2;        /* air potential temperature (k) at height zlvl above ground */
    double          sfctmp;
} estate_struct;

typedef struct eflux_struct
{
    double          solnet;     /* net downward solar radiation ((w m-2; positive) */
    double          etp;        /* potential evaporation (w m-2) */
    double          epsca;      /* */
    double          ssoil;      /* soil heat flux (w m-2: negative if downward from surface) */
    double          eta;        /* actual latent heat flux (w m-2: negative, if up from surface) */
    double          sheat;      /* sensible heat flux (w m-2: negative, if upward from surface) */
    double          fdown;      /* radiation forcing at the surface (w m-2) = soldn*(1-alb)+lwdn */
    double          lwdn;
    double          ec;         /* canopy water evaporation (w m-2) */
    double          edir;       /* direct soil evaporation (w m-2) */
    double          et[MAXLYR]; /* plant transpiration from a particular root (soil) layer (w m-2) */
    double          ett;        /* total plant transpiration (w m-2) */
    double          esnow;      /* sublimation from (or deposition to if <0) snowpack (w m-2) */
    double          soldn;      /* solar downward radiation (w m-2; positive, not net solar) */
    double          longwave;
    double          flx1;       /* precip-snow sfc (w m-2) */
    double          flx2;       /* freezing rain latent heat flux (w m-2) */
    double          flx3;       /* phase-change heat flux from snowmelt (w m-2) */
    double          solardirect;        /* direct component of downward solar radiation (w m-2) (not used) */
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
#endif

typedef struct elemic_struct
{
    double          intcp;
    double          sneqv;
    double          surf;
    double          unsat;
    double          gw;
#ifdef _NOAH_
    double          t1;
    double          snowh;
    double          stc[MAXLYR];
    double          smc[MAXLYR];
    double          sh2o[MAXLYR];
#endif
} elemic_struct;

#ifdef _BGC_
typedef struct restart_data_struct
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
    double          sminn;
    double          retransn;
    double          npool;
    double          day_leafc_litfall_increment;
    double          day_frootc_litfall_increment;
    double          day_livestemc_turnover_increment;
    double          day_livecrootc_turnover_increment;
    double          annmax_leafc;
    double          annmax_frootc;
    double          annmax_livestemc;
    double          annmax_livecrootc;
    double          dsr;
    double          dormant_flag;
    double          onset_flag;
    double          onset_counter;
    double          onset_gddflag;
    double          onset_fdd;
    double          onset_gdd;
    double          onset_swi;
    double          offset_flag;
    double          offset_counter;
    double          offset_fdd;
    double          offset_swi;
} restart_data_struct;
#endif

#ifdef _CYCLES_
typedef struct weather_struct
{
    //double          siteAltitude;
    //double          siteLatitude;
    //double          screeningHeight;
    //int             length;
    //double         *yearlyAmplitude;
    //double         *annualAverageTemperature;
    int            *lastDoy;
    double        **wind;
    //double        **ETref;
    //double        **precipitation;
    //double        **RHmax;
    //double        **RHmin;
    double        **solarRadiation;
    double        **tMax;
    double        **tMin;
    double        **vpd;
    double          atmosphericPressure;
} weather_struct;

typedef struct snow_struct
{
    double              snowCover;
} snow_struct;

typedef struct solute_struct
{
    double              soluteMass[MAXLYR];
    double              soluteMassAdsorbed[MAXLYR];
    double              soluteConc[MAXLYR];
    double              soluteFluxLat[4][MAXLYR];
    double              soluteFluxVert[MAXLYR];
} solute_struct;
#endif

#ifdef _BGC_
/* a structure to hold information on the annual co2 concentration */
typedef struct
{
    int             varco2;     /* (flag) 0=const 1=use file 2=const,file for Ndep */
    double          co2ppm;     /* (ppm)  constant CO2 concentration */
    double         *co2ppm_array;       /* (ppm)  annual CO2 concentration array */
    int            *co2year_array;      /* (year) year corresponding to the concentration value in co2ppm_arry */
    int             co2vals;    /* (num)  The number of CO2 concentration values in the co2ppm_array */
} co2control_struct;

/* a structure to hold annual nitrogen deposition data */
typedef struct
{
    int             varndep;    /* (flag) 0=const 1=use file  */
    double         *ndep_array; /* (kgN m-2 yr-1)  annual ndep array */
    int            *ndepyear_array;     /* (year) year corresponding to the ndep value in ndep_array */
    int             ndepvals;   /* (num)  The number of ndep values in the ndep_array */
    double          ndep;       /* (kgN/m2/yr) wet+dry atmospheric deposition of N */
    double          nfix;       /* (kgN/m2/yr) symbiotic+asymbiotic fixation of N */
} ndepcontrol_struct;

/* meteorological variable arrays */

/* inputs from mtclim, except for tavg and tavg_ra
 * which are used for an 11-day running average of daily average air T,
 * computed for the whole length of the met array prior to the 
 * daily model loop */
typedef struct
{
    double         *tmax;       /* (deg C) daily maximum air temperature */
    double         *tmin;       /* (deg C) daily minimum air temperature */
    double         *prcp;       /* (cm)    precipitation */
    double         *vpd;        /* (Pa)    vapor pressure deficit */
    double         *q2d;        /* (m3/m3) mixing ratio deficit */
    double         *swavgfd;    /* (W/m2)  daylight avg shortwave flux density */
    double         *par;        /* (W/m2)  photosynthetically active radiation */
    double         *dayl;       /* (s)     daylength */
    double         *prev_dayl;
    double         *tavg;       /* (deg C) daily average temperature */
    double         *tday;
    double         *tnight;
    double         *tsoil;
    double         *swc;
    double         *pa;
    double         *tavg_ra;    /* (deg C) 11-day running avg of daily avg temp */
    double         *latflux[4];
    double         *soilw;
    double         *sw_alb;
    double         *gl_bl;
    int            *flag;
} metarr_struct;

/* daily values that are passed to daily model subroutines */
typedef struct
{
    double          prcp;       /* (kg/m2) precipitation */
    double          tmax;       /* (deg C) daily maximum air temperature */
    double          tmin;       /* (deg C) daily minimum air temperature */
    double          tavg;       /* (deg C) daily average air temperature */
    double          tday;       /* (deg C) daylight average air temperature */
    double          tnight;     /* (deg C) nightime average air temperature */
    double          tsoil;      /* (deg C) daily soil temperature, avg, top 10 cm */
    double          swc;
    double          soilw;
    double          latflux[4];
    double          sw_alb;
    double          gl_bl;
    double          vpd;        /* (Pa)    vapor pressure deficit */
    double          q2d;        /* (m3/m3) mixing ratio deficit */
    double          swavgfd;    /* (W/m2)  daylight average shortwave flux */
    double          swabs;      /* (W/m2)  canopy absorbed shortwave flux */
    double          swtrans;    /* (W/m2)  transmitted shortwave flux */
    double          swabs_per_plaisun;  /* (W/m2) swabs per unit sunlit proj LAI */
    double          swabs_per_plaishade;        /* (W/m2) swabs per unit shaded proj LAI */
    double          ppfd_per_plaisun;   /* (umol/m2/s) ppfd per unit sunlit proj LAI */
    double          ppfd_per_plaishade; /* (umol/m2/s) ppfd per unit shaded proj LAI */
    double          par;        /* (W/m2)  photosynthetically active radiation */
    double          parabs;     /* (W/m2)  PAR absorbed by canopy */
    double          pa;         /* (Pa)    atmospheric pressure */
    double          co2;        /* (ppm)   atmospheric concentration of CO2 */
    double          dayl;       /* (s)     daylength */
    double          prev_dayl;  /* daylength from previous timestep (seconds) */
} metvar_struct;

/* carbon state initialization structure */
typedef struct
{
    double          max_leafc;  /* (kgC/m2) first-year displayed + stored leafc */
    double          max_stemc;  /* (kgC/m2) first-year total stem carbon */
} cinit_struct;

/* carbon state variables (including sums for sources and sinks) */
typedef struct
{
    double          leafc;      /* (kgC/m2) leaf C */
    double          leafc_storage;      /* (kgC/m2) leaf C storage */
    double          leafc_transfer;     /* (kgC/m2) leaf C transfer */
    double          frootc;     /* (kgC/m2) fine root C */
    double          frootc_storage;     /* (kgC/m2) fine root C storage */
    double          frootc_transfer;    /* (kgC/m2) fine root C transfer */
    double          livestemc;  /* (kgC/m2) live stem C */
    double          livestemc_storage;  /* (kgC/m2) live stem C storage */
    double          livestemc_transfer; /* (kgC/m2) live stem C transfer */
    double          deadstemc;  /* (kgC/m2) dead stem C */
    double          deadstemc_storage;  /* (kgC/m2) dead stem C storage */
    double          deadstemc_transfer; /* (kgC/m2) dead stem C transfer */
    double          livecrootc; /* (kgC/m2) live coarse root C */
    double          livecrootc_storage; /* (kgC/m2) live coarse root C storage */
    double          livecrootc_transfer;        /* (kgC/m2) live coarse root C transfer */
    double          deadcrootc; /* (kgC/m2) dead coarse root C */
    double          deadcrootc_storage; /* (kgC/m2) dead coarse root C storage */
    double          deadcrootc_transfer;        /* (kgC/m2) dead coarse root C transfer */
    double          gresp_storage;      /* (kgC/m2) growth respiration storage */
    double          gresp_transfer;     /* (kgC/m2) growth respiration transfer */
    double          cwdc;       /* (kgC/m2) coarse woody debris C */
    double          litr1c;     /* (kgC/m2) litter labile C */
    double          litr2c;     /* (kgC/m2) litter unshielded cellulose C */
    double          litr3c;     /* (kgC/m2) litter shielded cellulose C */
    double          litr4c;     /* (kgC/m2) litter lignin C */
    double          soil1c;     /* (kgC/m2) microbial recycling pool C (fast) */
    double          soil2c;     /* (kgC/m2) microbial recycling pool C (medium) */
    double          soil3c;     /* (kgC/m2) microbial recycling pool C (slow) */
    double          soil4c;     /* (kgC/m2) recalcitrant SOM C (humus, slowest) */
    double          cpool;      /* (kgC/m2) temporary photosynthate C pool */
    double          psnsun_src; /* (kgC/m2) SUM of gross PSN from sulit canopy */
    double          psnshade_src;       /* (kgC/m2) SUM of gross PSN from shaded canopy */
    double          leaf_mr_snk;        /* (kgC/m2) SUM of leaf maint resp */
    double          leaf_gr_snk;        /* (kgC/m2) SUM of leaf growth resp */
    double          froot_mr_snk;       /* (kgC/m2) SUM of fine root maint resp */
    double          froot_gr_snk;       /* (kgC/m2) SUM of fine root growth resp */
    double          livestem_mr_snk;    /* (kgC/m2) SUM of live stem maint resp */
    double          livestem_gr_snk;    /* (kgC/m2) SUM of live stem growth resp */
    double          deadstem_gr_snk;    /* (kgC/m2) SUM of dead stem growth resp */
    double          livecroot_mr_snk;   /* (kgC/m2) SUM of live coarse root maint resp */
    double          livecroot_gr_snk;   /* (kgC/m2) SUM of live coarse root growth resp */
    double          deadcroot_gr_snk;   /* (kgC/m2) SUM of dead coarse root growth resp */
    double          litr1_hr_snk;       /* (kgC/m2) SUM of labile litr microbial resp */
    double          litr2_hr_snk;       /* (kgC/m2) SUM of cellulose litr microbial resp */
    double          litr4_hr_snk;       /* (kgC/m2) SUM of lignin litr microbial resp */
    double          soil1_hr_snk;       /* (kgC/m2) SUM of fast microbial respiration */
    double          soil2_hr_snk;       /* (kgC/m2) SUM of medium microbial respiration */
    double          soil3_hr_snk;       /* (kgC/m2) SUM of slow microbial respiration */
    double          soil4_hr_snk;       /* (kgC/m2) SUM of recalcitrant SOM respiration */
    double          fire_snk;   /* (kgC/m2) SUM of fire losses */
} cstate_struct;

/* daily carbon flux variables */
typedef struct
{
    /* mortality fluxes */
    double          m_leafc_to_litr1c;  /* (kgC/m2/d) */
    double          m_leafc_to_litr2c;  /* (kgC/m2/d) */
    double          m_leafc_to_litr3c;  /* (kgC/m2/d) */
    double          m_leafc_to_litr4c;  /* (kgC/m2/d) */
    double          m_frootc_to_litr1c; /* (kgC/m2/d) */
    double          m_frootc_to_litr2c; /* (kgC/m2/d) */
    double          m_frootc_to_litr3c; /* (kgC/m2/d) */
    double          m_frootc_to_litr4c; /* (kgC/m2/d) */
    double          m_leafc_storage_to_litr1c;  /* (kgC/m2/d) */
    double          m_frootc_storage_to_litr1c; /* (kgC/m2/d) */
    double          m_livestemc_storage_to_litr1c;      /* (kgC/m2/d) */
    double          m_deadstemc_storage_to_litr1c;      /* (kgC/m2/d) */
    double          m_livecrootc_storage_to_litr1c;     /* (kgC/m2/d) */
    double          m_deadcrootc_storage_to_litr1c;     /* (kgC/m2/d) */
    double          m_leafc_transfer_to_litr1c; /* (kgC/m2/d) */
    double          m_frootc_transfer_to_litr1c;        /* (kgC/m2/d) */
    double          m_livestemc_transfer_to_litr1c;     /* (kgC/m2/d) */
    double          m_deadstemc_transfer_to_litr1c;     /* (kgC/m2/d) */
    double          m_livecrootc_transfer_to_litr1c;    /* (kgC/m2/d) */
    double          m_deadcrootc_transfer_to_litr1c;    /* (kgC/m2/d) */
    double          m_livestemc_to_cwdc;        /* (kgC/m2/d) */
    double          m_deadstemc_to_cwdc;        /* (kgC/m2/d) */
    double          m_livecrootc_to_cwdc;       /* (kgC/m2/d) */
    double          m_deadcrootc_to_cwdc;       /* (kgC/m2/d) */
    double          m_gresp_storage_to_litr1c;  /* (kgC/m2/d) */
    double          m_gresp_transfer_to_litr1c; /* (kgC/m2/d) */
    /* fire fluxes */
    double          m_leafc_to_fire;    /* (kgC/m2/d) */
    double          m_frootc_to_fire;   /* (kgC/m2/d) */
    double          m_leafc_storage_to_fire;    /* (kgC/m2/d) */
    double          m_frootc_storage_to_fire;   /* (kgC/m2/d) */
    double          m_livestemc_storage_to_fire;        /* (kgC/m2/d) */
    double          m_deadstemc_storage_to_fire;        /* (kgC/m2/d) */
    double          m_livecrootc_storage_to_fire;       /* (kgC/m2/d) */
    double          m_deadcrootc_storage_to_fire;       /* (kgC/m2/d) */
    double          m_leafc_transfer_to_fire;   /* (kgC/m2/d) */
    double          m_frootc_transfer_to_fire;  /* (kgC/m2/d) */
    double          m_livestemc_transfer_to_fire;       /* (kgC/m2/d) */
    double          m_deadstemc_transfer_to_fire;       /* (kgC/m2/d) */
    double          m_livecrootc_transfer_to_fire;      /* (kgC/m2/d) */
    double          m_deadcrootc_transfer_to_fire;      /* (kgC/m2/d) */
    double          m_livestemc_to_fire;        /* (kgC/m2/d) */
    double          m_deadstemc_to_fire;        /* (kgC/m2/d) */
    double          m_livecrootc_to_fire;       /* (kgC/m2/d) */
    double          m_deadcrootc_to_fire;       /* (kgC/m2/d) */
    double          m_gresp_storage_to_fire;    /* (kgC/m2/d) */
    double          m_gresp_transfer_to_fire;   /* (kgC/m2/d) */
    double          m_litr1c_to_fire;   /* (kgC/m2/d) */
    double          m_litr2c_to_fire;   /* (kgC/m2/d) */
    double          m_litr3c_to_fire;   /* (kgC/m2/d) */
    double          m_litr4c_to_fire;   /* (kgC/m2/d) */
    double          m_cwdc_to_fire;     /* (kgC/m2/d) */
    /* phenology fluxes from transfer pool */
    double          leafc_transfer_to_leafc;    /* (kgC/m2/d) */
    double          frootc_transfer_to_frootc;  /* (kgC/m2/d) */
    double          livestemc_transfer_to_livestemc;    /* (kgC/m2/d) */
    double          deadstemc_transfer_to_deadstemc;    /* (kgC/m2/d) */
    double          livecrootc_transfer_to_livecrootc;  /* (kgC/m2/d) */
    double          deadcrootc_transfer_to_deadcrootc;  /* (kgC/m2/d) */
    /* leaf and fine root litterfall */
    double          leafc_to_litr1c;    /* (kgC/m2/d) */
    double          leafc_to_litr2c;    /* (kgC/m2/d) */
    double          leafc_to_litr3c;    /* (kgC/m2/d) */
    double          leafc_to_litr4c;    /* (kgC/m2/d) */
    double          frootc_to_litr1c;   /* (kgC/m2/d) */
    double          frootc_to_litr2c;   /* (kgC/m2/d) */
    double          frootc_to_litr3c;   /* (kgC/m2/d) */
    double          frootc_to_litr4c;   /* (kgC/m2/d) */
    /* maintenance respiration fluxes */
    double          leaf_day_mr;        /* (kgC/m2/d) */
    double          leaf_night_mr;      /* (kgC/m2/d) */
    double          froot_mr;   /* (kgC/m2/d) */
    double          livestem_mr;        /* (kgC/m2/d) */
    double          livecroot_mr;       /* (kgC/m2/d) */
    /* photosynthesis flux */
    double          psnsun_to_cpool;    /* (kgC/m2/d) */
    double          psnshade_to_cpool;  /* (kgC/m2/d) */
    /* litter decomposition fluxes */
    double          cwdc_to_litr2c;     /* (kgC/m2/d) */
    double          cwdc_to_litr3c;     /* (kgC/m2/d) */
    double          cwdc_to_litr4c;     /* (kgC/m2/d) */
    double          litr1_hr;   /* (kgC/m2/d) */
    double          litr1c_to_soil1c;   /* (kgC/m2/d) */
    double          litr2_hr;   /* (kgC/m2/d) */
    double          litr2c_to_soil2c;   /* (kgC/m2/d) */
    double          litr3c_to_litr2c;   /* (kgC/m2/d) */
    double          litr4_hr;   /* (kgC/m2/d) */
    double          litr4c_to_soil3c;   /* (kgC/m2/d) */
    double          soil1_hr;   /* (kgC/m2/d) */
    double          soil1c_to_soil2c;   /* (kgC/m2/d) */
    double          soil2_hr;   /* (kgC/m2/d) */
    double          soil2c_to_soil3c;   /* (kgC/m2/d) */
    double          soil3_hr;   /* (kgC/m2/d) */
    double          soil3c_to_soil4c;   /* (kgC/m2/d) */
    double          soil4_hr;   /* (kgC/m2/d) */
    /* daily allocation fluxes from current GPP */
    double          cpool_to_leafc;     /* (kgC/m2/d) */
    double          cpool_to_leafc_storage;     /* (kgC/m2/d) */
    double          cpool_to_frootc;    /* (kgC/m2/d) */
    double          cpool_to_frootc_storage;    /* (kgC/m2/d) */
    double          cpool_to_livestemc; /* (kgC/m2/d) */
    double          cpool_to_livestemc_storage; /* (kgC/m2/d) */
    double          cpool_to_deadstemc; /* (kgC/m2/d) */
    double          cpool_to_deadstemc_storage; /* (kgC/m2/d) */
    double          cpool_to_livecrootc;        /* (kgC/m2/d) */
    double          cpool_to_livecrootc_storage;        /* (kgC/m2/d) */
    double          cpool_to_deadcrootc;        /* (kgC/m2/d) */
    double          cpool_to_deadcrootc_storage;        /* (kgC/m2/d) */
    double          cpool_to_gresp_storage;     /* (kgC/m2/d) */
    /* daily growth respiration fluxes */
    double          cpool_leaf_gr;      /* (kgC/m2/d) */
    double          cpool_leaf_storage_gr;      /* (kgC/m2/d) */
    double          transfer_leaf_gr;   /* (kgC/m2/d) */
    double          cpool_froot_gr;     /* (kgC/m2/d) */
    double          cpool_froot_storage_gr;     /* (kgC/m2/d) */
    double          transfer_froot_gr;  /* (kgC/m2/d) */
    double          cpool_livestem_gr;  /* (kgC/m2/d) */
    double          cpool_livestem_storage_gr;  /* (kgC/m2/d) */
    double          transfer_livestem_gr;       /* (kgC/m2/d) */
    double          cpool_deadstem_gr;  /* (kgC/m2/d) */
    double          cpool_deadstem_storage_gr;  /* (kgC/m2/d) */
    double          transfer_deadstem_gr;       /* (kgC/m2/d) */
    double          cpool_livecroot_gr; /* (kgC/m2/d) */
    double          cpool_livecroot_storage_gr; /* (kgC/m2/d) */
    double          transfer_livecroot_gr;      /* (kgC/m2/d) */
    double          cpool_deadcroot_gr; /* (kgC/m2/d) */
    double          cpool_deadcroot_storage_gr; /* (kgC/m2/d) */
    double          transfer_deadcroot_gr;      /* (kgC/m2/d) */
    /* annual turnover of storage to transfer pools */
    double          leafc_storage_to_leafc_transfer;    /* (kgC/m2/d) */
    double          frootc_storage_to_frootc_transfer;  /* (kgC/m2/d) */
    double          livestemc_storage_to_livestemc_transfer;    /* (kgC/m2/d) */
    double          deadstemc_storage_to_deadstemc_transfer;    /* (kgC/m2/d) */
    double          livecrootc_storage_to_livecrootc_transfer;  /* (kgC/m2/d) */
    double          deadcrootc_storage_to_deadcrootc_transfer;  /* (kgC/m2/d) */
    double          gresp_storage_to_gresp_transfer;    /* (kgC/m2/d) */
    /* turnover of live wood to dead wood */
    double          livestemc_to_deadstemc;     /* (kgC/m2/d) */
    double          livecrootc_to_deadcrootc;   /* (kgC/m2/d) */
} cflux_struct;

/* nitrogen state variables (including sums for sources and sinks) */
typedef struct
{
    double          leafn;      /* (kgN/m2) leaf N */
    double          leafn_storage;      /* (kgN/m2) leaf N */
    double          leafn_transfer;     /* (kgN/m2) leaf N */
    double          frootn;     /* (kgN/m2) fine root N */
    double          frootn_storage;     /* (kgN/m2) fine root N */
    double          frootn_transfer;    /* (kgN/m2) fine root N */
    double          livestemn;  /* (kgN/m2) live stem N */
    double          livestemn_storage;  /* (kgN/m2) live stem N */
    double          livestemn_transfer; /* (kgN/m2) live stem N */
    double          deadstemn;  /* (kgN/m2) dead stem N */
    double          deadstemn_storage;  /* (kgN/m2) dead stem N */
    double          deadstemn_transfer; /* (kgN/m2) dead stem N */
    double          livecrootn; /* (kgN/m2) live coarse root N */
    double          livecrootn_storage; /* (kgN/m2) live coarse root N */
    double          livecrootn_transfer;        /* (kgN/m2) live coarse root N */
    double          deadcrootn; /* (kgN/m2) dead coarse root N */
    double          deadcrootn_storage; /* (kgN/m2) dead coarse root N */
    double          deadcrootn_transfer;        /* (kgN/m2) dead coarse root N */
    double          cwdn;       /* (kgN/m2) coarse woody debris N */
    double          litr1n;     /* (kgN/m2) litter labile N */
    double          litr2n;     /* (kgN/m2) litter unshielded cellulose N */
    double          litr3n;     /* (kgN/m2) litter shielded cellulose N */
    double          litr4n;     /* (kgN/m2) litter lignin N */
    double          soil1n;     /* (kgN/m2) microbial recycling pool N (fast) */
    double          soil2n;     /* (kgN/m2) microbial recycling pool N (medium) */
    double          soil3n;     /* (kgN/m2) microbial recycling pool N (slow) */
    double          soil4n;     /* (kgN/m2) recalcitrant SOM N (humus, slowest) */
    double          sminn;      /* (kgN/m2) soil mineral N */
    double          retransn;   /* (kgN/m2) plant pool of retranslocated N */
    double          npool;      /* (kgN/m2) temporary plant N pool */
    double          nfix_src;   /* (kgN/m2) SUM of biological N fixation */
    double          ndep_src;   /* (kgN/m2) SUM of N deposition inputs */
    double          nleached_snk;       /* (kgN/m2) SUM of N leached */
    double          nvol_snk;   /* (kgN/m2) SUM of N lost to volatilization */
    double          fire_snk;   /* (kgN/m2) SUM of N lost to fire */
} nstate_struct;

/* daily nitrogen flux variables */
typedef struct
{
    /* mortality fluxes */
    double          m_leafn_to_litr1n;  /* (kgN/m2/d) */
    double          m_leafn_to_litr2n;  /* (kgN/m2/d) */
    double          m_leafn_to_litr3n;  /* (kgN/m2/d) */
    double          m_leafn_to_litr4n;  /* (kgN/m2/d) */
    double          m_frootn_to_litr1n; /* (kgN/m2/d) */
    double          m_frootn_to_litr2n; /* (kgN/m2/d) */
    double          m_frootn_to_litr3n; /* (kgN/m2/d) */
    double          m_frootn_to_litr4n; /* (kgN/m2/d) */
    double          m_leafn_storage_to_litr1n;  /* (kgN/m2/d) */
    double          m_frootn_storage_to_litr1n; /* (kgN/m2/d) */
    double          m_livestemn_storage_to_litr1n;      /* (kgN/m2/d) */
    double          m_deadstemn_storage_to_litr1n;      /* (kgN/m2/d) */
    double          m_livecrootn_storage_to_litr1n;     /* (kgN/m2/d) */
    double          m_deadcrootn_storage_to_litr1n;     /* (kgN/m2/d) */
    double          m_leafn_transfer_to_litr1n; /* (kgN/m2/d) */
    double          m_frootn_transfer_to_litr1n;        /* (kgN/m2/d) */
    double          m_livestemn_transfer_to_litr1n;     /* (kgN/m2/d) */
    double          m_deadstemn_transfer_to_litr1n;     /* (kgN/m2/d) */
    double          m_livecrootn_transfer_to_litr1n;    /* (kgN/m2/d) */
    double          m_deadcrootn_transfer_to_litr1n;    /* (kgN/m2/d) */
    double          m_livestemn_to_litr1n;      /* (kgN/m2/d) */
    double          m_livestemn_to_cwdn;        /* (kgN/m2/d) */
    double          m_deadstemn_to_cwdn;        /* (kgN/m2/d) */
    double          m_livecrootn_to_litr1n;     /* (kgN/m2/d) */
    double          m_livecrootn_to_cwdn;       /* (kgN/m2/d) */
    double          m_deadcrootn_to_cwdn;       /* (kgN/m2/d) */
    double          m_retransn_to_litr1n;       /* (kgN/m2/d) */
    /* fire fluxes */
    double          m_leafn_to_fire;    /* (kgN/m2/d) */
    double          m_frootn_to_fire;   /* (kgN/m2/d) */
    double          m_leafn_storage_to_fire;    /* (kgN/m2/d) */
    double          m_frootn_storage_to_fire;   /* (kgN/m2/d) */
    double          m_livestemn_storage_to_fire;        /* (kgN/m2/d) */
    double          m_deadstemn_storage_to_fire;        /* (kgN/m2/d) */
    double          m_livecrootn_storage_to_fire;       /* (kgN/m2/d) */
    double          m_deadcrootn_storage_to_fire;       /* (kgN/m2/d) */
    double          m_leafn_transfer_to_fire;   /* (kgN/m2/d) */
    double          m_frootn_transfer_to_fire;  /* (kgN/m2/d) */
    double          m_livestemn_transfer_to_fire;       /* (kgN/m2/d) */
    double          m_deadstemn_transfer_to_fire;       /* (kgN/m2/d) */
    double          m_livecrootn_transfer_to_fire;      /* (kgN/m2/d) */
    double          m_deadcrootn_transfer_to_fire;      /* (kgN/m2/d) */
    double          m_livestemn_to_fire;        /* (kgN/m2/d) */
    double          m_deadstemn_to_fire;        /* (kgN/m2/d) */
    double          m_livecrootn_to_fire;       /* (kgN/m2/d) */
    double          m_deadcrootn_to_fire;       /* (kgN/m2/d) */
    double          m_retransn_to_fire; /* (kgN/m2/d) */
    double          m_litr1n_to_fire;   /* (kgN/m2/d) */
    double          m_litr2n_to_fire;   /* (kgN/m2/d) */
    double          m_litr3n_to_fire;   /* (kgN/m2/d) */
    double          m_litr4n_to_fire;   /* (kgN/m2/d) */
    double          m_cwdn_to_fire;     /* (kgN/m2/d) */
    /* phenology fluxes from transfer pool */
    double          leafn_transfer_to_leafn;    /* (kgN/m2/d) */
    double          frootn_transfer_to_frootn;  /* (kgN/m2/d) */
    double          livestemn_transfer_to_livestemn;    /* (kgN/m2/d) */
    double          deadstemn_transfer_to_deadstemn;    /* (kgN/m2/d) */
    double          livecrootn_transfer_to_livecrootn;  /* (kgN/m2/d) */
    double          deadcrootn_transfer_to_deadcrootn;  /* (kgN/m2/d) */
    /* litterfall fluxes */
    double          leafn_to_litr1n;    /* (kgN/m2/d) */
    double          leafn_to_litr2n;    /* (kgN/m2/d) */
    double          leafn_to_litr3n;    /* (kgN/m2/d) */
    double          leafn_to_litr4n;    /* (kgN/m2/d) */
    double          leafn_to_retransn;  /* (kgN/m2/d) */
    double          frootn_to_litr1n;   /* (kgN/m2/d) */
    double          frootn_to_litr2n;   /* (kgN/m2/d) */
    double          frootn_to_litr3n;   /* (kgN/m2/d) */
    double          frootn_to_litr4n;   /* (kgN/m2/d) */
    /* deposition flux */
    double          ndep_to_sminn;      /* (kgN/m2/d) */
    double          nfix_to_sminn;      /* (kgN/m2/d) */
    /* litter and soil decomposition fluxes */
    double          cwdn_to_litr2n;     /* (kgN/m2/d) */
    double          cwdn_to_litr3n;     /* (kgN/m2/d) */
    double          cwdn_to_litr4n;     /* (kgN/m2/d) */
    double          litr1n_to_soil1n;   /* (kgN/m2/d) */
    double          sminn_to_soil1n_l1; /* (kgN/m2/d) */
    double          litr2n_to_soil2n;   /* (kgN/m2/d) */
    double          sminn_to_soil2n_l2; /* (kgN/m2/d) */
    double          litr3n_to_litr2n;   /* (kgN/m2/d) */
    double          litr4n_to_soil3n;   /* (kgN/m2/d) */
    double          sminn_to_soil3n_l4; /* (kgN/m2/d) */
    double          soil1n_to_soil2n;   /* (kgN/m2/d) */
    double          sminn_to_soil2n_s1; /* (kgN/m2/d) */
    double          soil2n_to_soil3n;   /* (kgN/m2/d) */
    double          sminn_to_soil3n_s2; /* (kgN/m2/d) */
    double          soil3n_to_soil4n;   /* (kgN/m2/d) */
    double          sminn_to_soil4n_s3; /* (kgN/m2/d) */
    double          soil4n_to_sminn;    /* (kgN/m2/d) */
    /* denitrification (volatilization) fluxes */
    double          sminn_to_nvol_l1s1; /* (kgN/m2/d) */
    double          sminn_to_nvol_l2s2; /* (kgN/m2/d) */
    double          sminn_to_nvol_l4s3; /* (kgN/m2/d) */
    double          sminn_to_nvol_s1s2; /* (kgN/m2/d) */
    double          sminn_to_nvol_s2s3; /* (kgN/m2/d) */
    double          sminn_to_nvol_s3s4; /* (kgN/m2/d) */
    double          sminn_to_nvol_s4;   /* (kgN/m2/d) */
    double          sminn_to_denitrif;  /* (kgN/m2/d) */

    /* leaching flux */
    double          sminn_leached;      /* (kgN/m2/d) */
    /* daily allocation fluxes */
    double          retransn_to_npool;  /* (kgN/m2/d) */
    double          sminn_to_npool;     /* (kgN/m2/d) */
    double          npool_to_leafn;     /* (kgN/m2/d) */
    double          npool_to_leafn_storage;     /* (kgN/m2/d) */
    double          npool_to_frootn;    /* (kgN/m2/d) */
    double          npool_to_frootn_storage;    /* (kgN/m2/d) */
    double          npool_to_livestemn; /* (kgN/m2/d) */
    double          npool_to_livestemn_storage; /* (kgN/m2/d) */
    double          npool_to_deadstemn; /* (kgN/m2/d) */
    double          npool_to_deadstemn_storage; /* (kgN/m2/d) */
    double          npool_to_livecrootn;        /* (kgN/m2/d) */
    double          npool_to_livecrootn_storage;        /* (kgN/m2/d) */
    double          npool_to_deadcrootn;        /* (kgN/m2/d) */
    double          npool_to_deadcrootn_storage;        /* (kgN/m2/d) */
    /* annual turnover of storage to transfer */
    double          leafn_storage_to_leafn_transfer;    /* (kgN/m2/d) */
    double          frootn_storage_to_frootn_transfer;  /* (kgN/m2/d) */
    double          livestemn_storage_to_livestemn_transfer;    /* (kgN/m2/d) */
    double          deadstemn_storage_to_deadstemn_transfer;    /* (kgN/m2/d) */
    double          livecrootn_storage_to_livecrootn_transfer;  /* (kgN/m2/d) */
    double          deadcrootn_storage_to_deadcrootn_transfer;  /* (kgN/m2/d) */
    /* turnover of live wood to dead wood, with retranslocation */
    double          livestemn_to_deadstemn;     /* (kgN/m2/d) */
    double          livestemn_to_retransn;      /* (kgN/m2/d) */
    double          livecrootn_to_deadcrootn;   /* (kgN/m2/d) */
    double          livecrootn_to_retransn;     /* (kgN/m2/d) */
} nflux_struct;

/* temporary nitrogen variables for reconciliation of decomposition
 * immobilization fluxes and plant growth N demands */
typedef struct
{
    double          mineralized;
    double          potential_immob;
    double          plitr1c_loss;
    double          pmnf_l1s1;
    double          plitr2c_loss;
    double          pmnf_l2s2;
    double          plitr4c_loss;
    double          pmnf_l4s3;
    double          psoil1c_loss;
    double          pmnf_s1s2;
    double          psoil2c_loss;
    double          pmnf_s2s3;
    double          psoil3c_loss;
    double          pmnf_s3s4;
    double          psoil4c_loss;
    double          kl4;
} ntemp_struct;

/* phenological control arrays */
typedef struct
{
    int            *remdays_curgrowth;  /* (nmetdays) days left in current growth season */
    int            *remdays_transfer;   /* (nmetdays) number of transfer days remaining */
    int            *remdays_litfall;    /* (nmetdays) number of litfall days remaining */
    int            *predays_transfer;   /* (nmetdays) number of transfer days previous */
    int            *predays_litfall;    /* (nmetdays) number of litfall days previous */
} phenarray_struct;

/* daily phenological data array */
typedef struct
{
    double          remdays_curgrowth;  /* days left in current growth season */
    double          remdays_transfer;   /* number of transfer days remaining */
    double          remdays_litfall;    /* number of litfall days remaining */
    double          predays_transfer;   /* number of transfer days previous */
    double          predays_litfall;    /* number of litfall days previous */
} phenology_struct;

/* ecophysiological variables */
typedef struct
{
    double          day_leafc_litfall_increment;        /* (kgC/m2/d) rate leaf litfall */
    double          day_frootc_litfall_increment;       /* (kgC/m2/d) rate froot litfall */
    double          day_livestemc_turnover_increment;   /* (kgC/m2/d) rate livestem turnover */
    double          day_livecrootc_turnover_increment;  /* (kgC/m2/d) rate livecroot turnover */
    double          annmax_leafc;       /* (kgC/m2) annual maximum daily leaf C */
    double          annmax_frootc;      /* (kgC/m2) annual maximum daily froot C */
    double          annmax_livestemc;   /* (kgC/m2) annual maximum daily livestem C */
    double          annmax_livecrootc;  /* (kgC/m2) annual maximum daily livecroot C */
    double          dsr;        /* (days) number of days since rain, for soil evap */
    double          proj_lai;   /* (DIM) live projected leaf area index */
    double          all_lai;    /* (DIM) live all-sided leaf area index */
    double          plaisun;    /* (DIM) sunlit projected leaf area index */
    double          plaishade;  /* (DIM) shaded projected leaf area index */
    double          sun_proj_sla;       /* (m2/kgC) sunlit projected SLA */
    double          shade_proj_sla;     /* (m2/kgC) shaded projected SLA */
    double          psi;        /* (MPa) water potential of soil and leaves */
    double          vwc;        /* (DIM) volumetric water content */
    double          dlmr_area_sun;      /* (umolC/m2projected leaf area/s) sunlit leaf MR */
    double          dlmr_area_shade;    /* (umolC/m2projected leaf area/s) shaded leaf MR */
    double          gl_t_wv_sun;        /* (m/s) leaf-scale conductance to transpired water */
    double          gl_t_wv_shade;      /* (m/s) leaf-scale conductance to transpired water */
    double          assim_sun;  /* (umol/m2/s) sunlit assimilation per unit pLAI */
    double          assim_shade;        /* (umol/m2/s) shaded assimilation per unit pLAI */
    /* decomp variables */
    double          t_scalar;   /* (DIM) decomp temperature scalar */
    double          w_scalar;   /* (DIM) decomp water scalar */
    double          rate_scalar;        /* (DIM) decomp combined scalar */
    double          daily_gross_nmin;   /* (kgN/m2/d) daily gross N mineralization */
    double          daily_gross_nimmob; /* (kgN/m2/d) daily gross N immobilization */
    double          daily_net_nmin;     /* (kgN/m2/d) daily net N mineralization */
    double          fpi;        /* (DIM) fraction of potential immobilization */

    /* the following are optional outputs, usually set if the appropriate
     * functions are called with the flag verbose = 1 */
    double          m_tmin;     /* (DIM) freezing night temperature multiplier */
    double          m_psi;      /* (DIM) water potential multiplier */
    double          m_co2;      /* (DIM) atmospheric [CO2] multiplier */
    double          m_ppfd_sun; /* (DIM) PAR flux density multiplier */
    double          m_ppfd_shade;       /* (DIM) PAR flux density multiplier */
    double          m_vpd;      /* (DIM) vapor pressure deficit multiplier */
    double          m_final_sun;        /* (DIM) product of all other multipliers */
    double          m_final_shade;      /* (DIM) product of all other multipliers */
    double          gl_bl;      /* (m/s) leaf boundary layer conductance */
    double          gl_c;       /* (m/s) leaf cuticular conductance */
    double          gl_s_sun;   /* (m/s) leaf-scale stomatal conductance */
    double          gl_s_shade; /* (m/s) leaf-scale stomatal conductance */
    double          gl_e_wv;    /* (m/s) leaf conductance to evaporated water */
    double          gl_sh;      /* (m/s) leaf conductance to sensible heat */
    double          gc_e_wv;    /* (m/s) canopy conductance to evaporated water */
    double          gc_sh;      /* (m/s) canopy conductance to sensible heat */

    /* diagnostic variables for ouput purposes only */
    double          ytd_maxplai;        /* (DIM) year-to-date maximum projected LAI */

    /* Variables below are added for PIHM-BGC */
    double          dormant_flag;       /* dormancy flag */
    double          days_active;        /* number of days since last dormancy */
    double          onset_flag; /* onset flag */
    double          onset_counter;      /* onset days counter */
    double          onset_gddflag;      /* onset flag for growing degree day sum */
    double          onset_fdd;  /* onset freezing degree days counter */
    double          onset_gdd;  /* onset growing degree days */
    double          onset_swi;  /* onset soil water index */
    double          offset_flag;        /* offset flag */
    double          offset_counter;     /* offset days counter */
    double          offset_fdd; /* offset freezing degree days counter */
    double          offset_swi; /* offset soil water index */
    double          lgsf;       /* long growing season factor [0-1] */
    double          bglfr;      /* background litterfall rate (1/s) */
    double          bgtr;       /* background transfer growth rate (1/s) */
    //   double dayl;     /* daylength (seconds) */
    double          annavg_t2m; /* annual average 2m air temperature (K) */
    double          tempavg_t2m;        /* temporary average 2m air temperature (K) */
    double          gpp;        /* GPP flux before downregulation (gC/m2/s) */
    double          availc;     /* C flux available for allocation (gC/m2/s) */
    double          xsmrpool_recover;   /* C flux assigned to recovery of negative cpool (gC/m2/s) */
    double          xsmrpool_c13ratio;  /* C13/C(12+13) ratio for xsmrpool (proportion) */
    double          alloc_pnow; /* fraction of current allocation to display as new growth (DIM) */
    double          c_allometry;        /* C allocation index (DIM) */
    double          n_allometry;        /* N allocation index (DIM) */
    double          plant_ndemand;      /* N flux required to support initial GPP (gN/m2/s) */
    double          tempsum_potential_gpp;      /* temporary annual sum of potential GPP */
    double          annsum_potential_gpp;       /* annual sum of potential GPP */
    double          tempmax_retransn;   /* temporary annual max of retranslocated N pool (gN/m2) */
    double          annmax_retransn;    /* annual max of retranslocated N pool (gN/m2) */
    double          avail_retransn;     /* N flux available from retranslocation pool (gN/m2/s) */
    double          plant_nalloc;       /* total allocated N flux (gN/m2/s) */
    double          plant_calloc;       /* total allocated C flux (gC/m2/s) */
    double          excess_cflux;       /* C flux not allocated due to downregulation (gC/m2/s) */
    double          downreg;    /* fractional reduction in GPP due to N limitation (DIM) */
    double          prev_leafc_to_litter;       /* previous timestep leaf C litterfall flux (gC/m2/s) */
    double          prev_frootc_to_litter;      /* previous timestep froot C litterfall flux (gC/m2/s) */
    double          tempsum_npp;        /* temporary annual sum of NPP (gC/m2/yr) */
    double          annsum_npp; /* annual sum of NPP (gC/m2/yr) */
    double          tempsum_litfall;    /* temporary annual sum of litfall (gC/m2/yr) */
    double          annsum_litfall;     /* annual sum of litfall (gC/m2/yr) */
    double          rc13_canair;        /* C13O2/C12O2 in canopy air */
    double          rc13_psnsun;        /* C13O2/C12O2 in sunlit canopy psn flux */
    double          rc13_psnsha;        /* C13O2/C12O2 in shaded canopy psn flux */

    double          old_c_balance;
    double          old_n_balance;
} epvar_struct;

typedef struct
{
    int             woody;      /* (flag) 1=woody, 0=non-woody */
    int             evergreen;  /* (flag) 1=evergreen, 0=deciduous */
    int             c3_flag;    /* (flag) 1 = C3,  0 = C4 */
    int             phenology_flag;     /* (flag) 1=phenology model, 0=user defined */
    int             onday;      /* (yday) yearday leaves on */
    int             offday;     /* (yday) yearday leaves off */
    double          transfer_days;      /* (prop.) fraction of growth period for transfer */
    double          litfall_days;       /* (prop.) fraction of growth period for litfall */
    double          leaf_turnover;      /* (1/yr) annual leaf turnover fraction */
    double          froot_turnover;     /* (1/yr) annual fine root turnover fraction */
    double          livewood_turnover;  /* (1/yr) annual live wood turnover fraction */
    double          daily_mortality_turnover;   /* (1/day) daily mortality turnover */
    double          daily_fire_turnover;        /* (1/day) daily fire turnover */
    double          alloc_frootc_leafc; /* (ratio) new fine root C to new leaf C */
    double          alloc_newstemc_newleafc;    /* (ratio) new stem C to new leaf C */
    double          alloc_newlivewoodc_newwoodc;        /* (ratio) new livewood C:new wood C */
    double          alloc_crootc_stemc; /* (ratio) new live croot C to new live stem C */
    double          alloc_prop_curgrowth;       /* (prop.) daily allocation to current growth */
    double          avg_proj_sla;       /* (m2/kgC) canopy average proj. SLA */
    double          sla_ratio;  /* (DIM) ratio of shaded to sunlit projected SLA */
    double          lai_ratio;  /* (DIM) ratio of (all-sided LA / one-sided LA) */
    double          int_coef;   /* (kg/kg/LAI/d) canopy precip interception coef */
    double          ext_coef;   /* (DIM) canopy light extinction coefficient */
    double          flnr;       /* (kg NRub/kg Nleaf) leaf N in Rubisco */
    double          psi_open;   /* (MPa) psi at start of conductance reduction */
    double          psi_close;  /* (MPa) psi at complete conductance reduction */
    double          vpd_open;   /* (Pa)  vpd at start of conductance reduction */
    double          vpd_close;  /* (Pa)  vpd at complete conductance reduction */
    double          gl_smax;    /* (m/s) maximum leaf-scale stomatal conductance */
    double          gl_c;       /* (m/s) leaf-scale cuticular conductance */
    double          gl_bl;      /* (m/s) leaf-scale boundary layer conductance */
    double          froot_cn;   /* (kgC/kgN) C:N for fine roots */
    double          leaf_cn;    /* (kgC/kgN) C:N for leaves */
    double          livewood_cn;        /* (kgC/kgN) C:N for live wood */
    double          deadwood_cn;        /* (kgC/kgN) C:N for dead wood */
    double          leaflitr_cn;        /* (kgC/kgN) constant C:N for leaf litter */
    double          leaflitr_flab;      /* (DIM) leaf litter labile fraction */
    double          leaflitr_fucel;     /* (DIM) leaf litter unshielded cellulose fract. */
    double          leaflitr_fscel;     /* (DIM) leaf litter shielded cellulose fract. */
    double          leaflitr_flig;      /* (DIM) leaf litter lignin fraction */
    double          frootlitr_flab;     /* (DIM) froot litter labile fraction */
    double          frootlitr_fucel;    /* (DIM) froot litter unshielded cellulose fract */
    double          frootlitr_fscel;    /* (DIM) froot litter shielded cellulose fract */
    double          frootlitr_flig;     /* (DIM) froot litter lignin fraction */
    double          deadwood_fucel;     /* (DIM) dead wood unshileded cellulose fraction */
    double          deadwood_fscel;     /* (DIM) dead wood shielded cellulose fraction */
    double          deadwood_flig;      /* (DIM) dead wood lignin fraction */

    //double          topt;
    //double          rgl;
    //double          hs;
    //double          smcref;
    //double          smcwlt;
} epconst_struct;

/* structure for the photosynthesis routine */
typedef struct
{
    int             c3;         /* (flag) set to 1 for C3 model, 0 for C4 model */
    double          pa;         /* (Pa) atmospheric pressure */
    double          co2;        /* (ppm) atmospheric [CO2] */
    double          t;          /* (deg C) temperature */
    double          lnc;        /* (kg Nleaf/m2) leaf N per unit sunlit leaf area */
    double          flnr;       /* (kg NRub/kg Nleaf) fract. of leaf N in Rubisco */
    double          ppfd;       /* (umol/m2/s) PAR flux per unit sunlit leaf area */
    double          g;          /* (umol/m2/s/Pa) conductance to CO2 */
    double          dlmr;       /* (umol/m2/s) day leaf m. resp, proj. area basis */
    double          Ci;         /* (Pa) intercellular [CO2] */
    double          O2;         /* (Pa) atmospheric [O2] */
    double          Ca;         /* (Pa) atmospheric [CO2] */
    double          gamma;      /* (Pa) CO2 compensation point, no Rd */
    double          Kc;         /* (Pa) MM constant carboxylation */
    double          Ko;         /* (Pa) MM constant oxygenation */
    double          Vmax;       /* (umol/m2/s) max rate carboxylation */
    double          Jmax;       /* (umol/m2/s) max rate electron transport */
    double          J;          /* (umol/m2/s) rate of RuBP regeneration */
    double          Av;         /* (umol/m2/s) carboxylation limited assimilation */
    double          Aj;         /* (umol/m2/s) RuBP regen limited assimilation */
    double          A;          /* (umol/m2/s) final assimilation rate */
} psn_struct;

typedef struct
{
    double          daily_npp;  /* kgC/m2/day = GPP - Rmaint - Rgrowth */
    double          daily_nep;  /* kgC/m2/day = NPP - Rheterotroph */
    double          daily_nee;  /* kgC/m2/day = NEP - fire losses */
    double          daily_gpp;  /* kgC/m2/day  gross PSN source */
    double          daily_mr;   /* kgC/m2/day  maintenance respiration */
    double          daily_gr;   /* kgC/m2/day  growth respiration */
    double          daily_hr;   /* kgC/m2/day  heterotrophic respiration */
    double          daily_fire; /* kgC/m2/day  fire losses */
    double          daily_litfallc;     /* kgC/m2/day  total litterfall */
    double          daily_et;   /* kgW/m2/day daily evapotranspiration */
    double          daily_evap; /* kgW/m2/day daily evaporation */
    double          daily_trans;        /* kgW/m2/day daily transpiration */
    double          daily_outflow;      /* kgW/m2/day daily outflow */
    double          daily_soilw;        /* kgW/m2/day daily soilw */
    double          daily_snoww;        /* kgW/m2/day daily snoww */
    double          cum_npp;    /* kgC/m2  Summed over entire simulation */
    double          cum_nep;    /* kgC/m2  Summed over entire simulation */
    double          cum_nee;    /* kgC/m2  Summed over entire simulation */
    double          cum_gpp;    /* kgC/m2  Summed over entire simulation */
    double          cum_mr;     /* kgC/m2  Summed over entire simulation */
    double          cum_gr;     /* kgC/m2  Summed over entire simulation */
    double          cum_hr;     /* kgC/m2  Summed over entire simulation */
    double          cum_fire;   /* kgC/m2  Summed over entire simulation */
    double          vegc;       /* kgC/m2  total vegetation C */
    double          litrc;      /* kgC/m2  total litter C */
    double          soilc;      /* kgC/m2  total soil C */
    double          totalc;     /* kgC/m2  total of vegc, litrc, and soilc */
} summary_struct;

typedef struct epclist_struct
{
    int             nvegtypes;  /* number of vegetation types */
    epconst_struct *epc;        /* pointer to array of epc structures */
} epclist_struct;
#endif

typedef struct elem_struct
{
    int             node[3];    /* Counterclock-wise */
    int             nabr[3];    /* neighbor i shares edge i
                                 * (0: on boundary) */
    int             ind;

    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
    elemforc_struct forc;
    elemic_struct   ic;

#ifdef _DAILY_
    daily_struct    daily;
#endif
    pstate_struct       ps;

    wstate_struct       ws;
    wstate_struct       ws0;
    wflux_struct       wf;
#ifdef _NOAH_
    wflux_struct       avgwf;

    estate_struct       es;
    eflux_struct       ef;
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
#endif
#ifdef _BGC_
    restart_data_struct restart_input;
    restart_data_struct restart_output;
    metarr_struct   metarr;     /* meteorological data array */
    metvar_struct   metv;
    cinit_struct    cinit;      /* first-year values for leafc and stemc */
    cstate_struct   cs;         /* carbon state variables */
    cflux_struct    cf;
    nstate_struct   ns;         /* nitrogen state variables */
    nflux_struct    nf;
    psn_struct      psn_sun;
    psn_struct      psn_shade;
    ntemp_struct    nt;
    summary_struct  summary;

    phenology_struct phen;

    epconst_struct  epc;        /* ecophysiological constants */
    epvar_struct    epv;
#endif

} elem_struct;

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

typedef struct riveric_struct
{
    double          stage;
    double          gw;
} riveric_struct;

typedef struct river_struct
{
    topo_struct     topo;
    shp_struct      shp;
    matl_struct     matl;
    elemforc_struct forc;
    wstate_struct       ws;
    wstate_struct       ws0;
    wflux_struct       wf;
    riveric_struct  ic;
#ifdef _DAILY_
    daily_struct    daily;
#endif
    int             leftele;    /* Left neighboring element */
    int             rightele;   /* Right neighboring element */
    int             fromnode;   /* Upstream Node no. */
    int             tonode;     /* Dnstream Node no. */
    int             down;       /* down stream segment */
#ifdef _CYCLES_
    solute_struct   NO3sol;
    solute_struct   NH4sol;
#endif
#ifdef _BGC_
    metarr_struct   metarr;     /* meteorological data array */
    double          sminn;
    double          nleached_snk;
    double          sminn_leached;
#endif
} river_struct;

typedef struct calib_struct
{
    double          ksath;
    double          ksatv;
    double          kinfv;
    double          kmach;
    double          kmacv;
    double          dinf;
    double          rzd;
    double          dmac;
    double          porosity;
    double          alpha;
    double          beta;
    double          areafv;
    double          areafh;
    double          vegfrac;
    double          albedo;
    double          rough;
    double          prcp;
    double          sfctmp;

    double          ec;
    double          ett;
    double          edir;

    double          rivrough;
    double          rivksath;
    double          rivksatv;
    double          rivbedthick;
    double          rivdepth;
    double          rivshpcoeff;

#ifdef _NOAH_
    double          thetaref;   /* ys */
    double          thetaw;     /* ys */
    double          rsmin;      /* ys */
    double          drip;       /* ys */
    double          intcp;      /* ys */
    double          czil;       /* ys */
    double          fxexp;      /* ys */
    double          cfactr;     /* ys */
    double          rgl;        /* ys */
    double          hs;         /* ys */
#endif
#ifdef _RT_
    double          pco2;
    double          keq;
    double          ssa;
    double          site_den;
    double          prep_conc;
#endif
} calib_struct;

typedef struct ctrl_struct
{
    int             ascii;      /* ys: add ascii output model
                                 * (default is binary */
    int             write_ic;   /* ys: runs model as spinup. model output at
                                 * the last step will be saved in .init */
    int             solver;     /* solver type */
    int             nstep;      /* number of external time steps (when
                                 * results can be printed) for the
                                 * whole simulation */
    int             nprint;     /* ys: number of variables for output */

    /* time interval to output average values of variables
     * variables will not be printed if time intervals are set to 0 */
    int             prtvrbl[NUM_PRINT];

    int             init_type;  /* initialization mode */

    int             unsat_mode; /* unsat mode */
    int             surf_mode;  /* surface overland flow mode */
    int             riv_mode;   /* river routing mode */


    double          abstol;     /* absolute tolerance */
    double          reltol;     /* relative tolerance */
    double          initstep;   /* initial step size */
    double          maxstep;    /* maximum step size */

    int             etstep;     /* step for et from interception */

    int             starttime;  /* start time of simulation */
    int             endtime;    /* end time of simulation */

    int             stepsize;

    int            *tout;
#ifdef _NOAH_
    int             nsoil;
    double          sldpth[MAXLYR];
    int             rad_mode;   /* radiation mode; 1: topographic, 0: uniform */
#endif

#ifdef _BGC_
    double          simstarttime;       /* start time of simulation */
    double          simendtime; /* end time of simulation */
    int             spinup;     /* (flag) 1=spinup run, 0=normal run */
    int             maxspinyears;       /* maximum number of years for spinup run */
    int             dodaily;    /* flag for daily output */
    int             domonavg;   /* flag for monthly average of daily outputs */
    int             doannavg;   /* flag for annual average of daily outputs */
    int             doannual;   /* flag for annual output */
    int             ndayout;    /* number of daily outputs */
    int             nannout;    /* number of annual outputs */
    int            *daycodes;   /* array of indices for daily outputs */
    int            *anncodes;   /* array of indices for annual outputs */
    int             read_restart;       /* flag to read restart file */
    int             write_restart;      /* flag to write restart file */
    int             keep_metyr; /* (flag) 1=retain restart metyr, 0=reset metyr */
    int             onscreen;   /* (flag) 1=show progress on-screen 0=don't */
    int             spinupstartyear;    /* first met year for spinup */
    int             spinupendyear;      /* last met year for spinup */
    int             spinupstart;        /* start time of spinup */
    int             spinupend;  /* end time of spinup */

    cstate_struct   cs;
    nstate_struct   ns;
    cinit_struct    cinit;
#endif
} ctrl_struct;

typedef struct prtctrl_struct
{
    char            name[MAXSTRING];
    int             intvl;
    int             nvrbl;
    double        **vrbl;
    double         *buffer;
} prtctrl_struct;

/* Model_data definition */
typedef struct pihm_struct
{
    int             numele;     /* number of elements */
    int             numriv;     /* number of rivere segments */

    double          longitude;
    double          latitude;
    double          elevation;

    /*
     * Input
     */
    /* File names */
    filename_struct filename;

    /* Input look-up talbes */
    meshtbl_struct  meshtbl;
    atttbl_struct   atttbl;
    soiltbl_struct  soiltbl;    /* Store Soil Information */
    geoltbl_struct  geoltbl;    /* Store Soil Information */
    lctbl_struct    lctbl;      /* Store Land Cover Information */
    rivtbl_struct   rivtbl;     /* Store River Segment Information */
    shptbl_struct   shptbl;     /* Store River Shape Information */
    matltbl_struct  matltbl;    /* Store River Bank Material Information */
#ifdef _NOAH_
    noahtbl_struct  noahtbl;
#endif
#ifdef _CYCLES_
    agtbl_struct    agtbl;
    croptbl_struct  croptbl;
    mgmttbl_struct  mgmttbl;
#endif
#ifdef _BGC_
    co2control_struct co2;      /* CO2 concentration information */
    ndepcontrol_struct ndepctrl;        /* Nitrogen deposition control structure */
    epclist_struct  epclist;
#endif
    forc_struct     forc;

    elem_struct    *elem;       /* Store Element Information */
    river_struct   *riv;

    calib_struct    cal;
    ctrl_struct     ctrl;
    prtctrl_struct  prtctrl[NUM_PRINT];
}              *pihm_struct;

#endif
