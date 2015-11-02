#ifndef CYCLES_STRUCT_HEADER
#define CYCLES_STRUCT_HEADER
/*****************************************************************************
 *
 * FILE NAME:   Cycles_struct.h
 * PURPOSE:     Defines Cycles data structure.
 *
 * Variables                Type        Description
 * ==========               ==========  ====================
 * ---------------------------------------------------------------------------
 * SimControlStruct         struct      Simulation control structure
 *  simStartYear            int         Simulation start year
 *  simEndYear              int         Simulation end year
 *  totalYears              int         Total simulation years
 *  yearsInRotation         int         Years in rotation
 *  adjustedYields          int
 *  hourlyInfiltration      int         Flag to indicate whether hourly
 *                                        infiltration is used
 *  automaticNitrogen       int         Flag to indicate whether automatic
 *                                        Nitrogen is used
 *  automaticPhosphorus     int             
 *  automaticSulfur         int         Flag to indicate
 *  cropDailyOutput         int         Flag to control daily crop output
 *  soilDailyOutput         int         Flag to control daily soil output
 *  nitrogenDailyOutput     int         Flag to control daily nitrogen output
 *  waterDailyOutput        int         Flag to control daily water output
 *  weatherDailyOutput      int         Flag to control daily weather output
 *  residueDailyOutput      int         Flag to control daily residue output
 *  soilCarbonDailyOutput   int         Flag to control daily soil carbon
 *                                        output
 *  annualSoilOutput        int         Flag to control annual soil output
 *  profileOutput           int         Flag to control profile output
 *  seasonOutput            int         Flag to control seasonal output
 * ---------------------------------------------------------------------------
 * SoilStruct               struct      Soil structure
 *  totalLayers             int         Total soil layers [input]
 *  Curve_Number            double      Curve number [input]
 *  Percent_Slope           double      Slope [input]
 *  annualTemperaturePhase  double      Annual soil temperature phase
 *  dampingDepth            double*     Soil temperature damping depth (m)
 *  cumulativeDepth         double*     Depth to the bottom of layer (m)
 *  nodeDepth               double*     Depth of node (m)
 *  layerThickness          double*     Measured layer thickness (m) [input]
 *  Clay                    double*     Clay fraction [input]
 *  Sand                    double*     Sand fraction [input]
 *  IOM                     double*     Initial Organic Matter [input]
 *  NO3                     double*     Initial Nitrate (kg/ha) [input]
 *  NH4                     double*     Initial Ammonium (kg/ha) [input]
 *  BD                      double*     Bulk Density (Mg/m3)
 *  FC                      double*     Field Capacity water content [input]
 *  PWP                     double*     Permanent Wilting Point [input]
 *  Porosity                double*     Saturation water content (m3/m3)
 *  PAW                     double*     Maximum plant available water
 *  FC_WaterPotential       double*     Estimate water potential at field
 *                                        capacity
 *  airEntryPotential       double*     Calculate Air Entry Potential
 *  B_Value                 double*     Calculated "B" value
 *  M_Value                 double*     Calculated "M" value
 *  n2o                     double*     Temporary output of n2o/layer
 *  SOC_Conc                double*     (g C/kg soil)
 *  SOC_Mass                double*     Soil organic carbon (Mg/ha)
 *  SON_Mass                double*     Soil organic Nitrogen (Mg/ha)
 *  MBC_Mass                double*     Microbial biomass C (Mg/ha)
 *  MBN_Mass                double*     Microbial Biomass N (Mg/ha)
 *  SOCProfile              double      Profile soil organic C
 *  SONProfile              double      Profile soil organic N
 *  C_Humified              double      Carbon humified from residues, roots,
 *                                        rizho, and manure
 *  C_ResidueRespired       double      Carbon respired from residues, roots,
 *                                        rizho, and manure
 *  C_SoilRespired          double      Carbon respired from soil organic
 *                                        carbon only
 *  soilTemperature         double*     Soil temperature (Degree Celsius)
 *  waterContent            double*     Volumetric water content (m3/m3)
 *  waterUptake             double*     Layer water uptake
 *  pH                      double*     Soil pH
 *  evaporationVol          double      Soil evaporation (mm)
 *  residueEvaporationVol   double      Evaporation from residue (mm)
 *  infiltrationVol         double      Soil infiltration (mm)
 *  runoffVol               double      Runoff (mm)
 *  irrigationVol           double      Irrigation (mm)
 *  drainageVol             double      Drainage (mm)
 *  NO3Leaching             double      NO3 leaching (Mg N/ha)
 *  NH4Leaching             double      NH4 leaching (Mg N/ha)
 *  NO3Profile              double
 *  NH4Profile              double
 *  N_Immobilization        double
 *  N_Mineralization        double
 *  N_NetMineralization     double
 *  NH4_Nitrification       double
 *  N2O_Nitrification       double
 *  NO3_Denitrification     double
 *  N2O_Denitrification     double
 *  NH4_Volatilization      double
 *
 *
 ***************************************************************************/
typedef struct SimControlStruct
{
    int             simStartYear;
    int             simEndYear;
    int             totalYears;
    int             yearsInRotation;

    int             adjustedYields;
    int             hourlyInfiltration;
    int             automaticNitrogen;
    int             automaticPhosphorus;
    int             automaticSulfur;
    int             cropDailyOutput;
    int             soilDailyOutput;
    int             nitrogenDailyOutput;
    int             waterDailyOutput;
    int             weatherDailyOutput;
    int             residueDailyOutput;
    int             soilCarbonDailyOutput;
    int             annualSoilOutput;
    int             profileOutput;
    int             seasonOutput;

    char            crop_filename[128];
    char            operation_filename[128];
    char            weather_filename[128];
    char            soil_filename[128];
} SimControlStruct;

typedef struct SoilStruct
{
    int             totalLayers;
    double          Curve_Number;
    double          Percent_Slope;

    double          annualTemperaturePhase;
    double          dampingDepth;

    double         *cumulativeDepth;
    double         *nodeDepth;
    double         *layerThickness;
    double         *Clay;
    double         *Sand;
    double         *IOM;
    double         *NO3;
    double         *NH4;
    double         *BD;
    double         *FC;
    double         *PWP;

    double         *Porosity;
    double         *PAW;
    double         *FC_WaterPotential;
    double         *airEntryPotential;
    double         *B_Value;
    double         *M_Value;

    double         *n2o;

    double         *SOC_Conc;
    double         *SOC_Mass;
    double         *SON_Mass;
    double         *MBC_Mass;
    double         *MBN_Mass;
    double          SOCProfile;
    double          SONProfile;

    double          C_Humified;
    double          C_ResidueRespired;
    double          C_SoilRespired;

    double         *soilTemperature;
    double         *waterContent;
    double         *waterUptake;
    double         *pH;

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
} SoilStruct;

typedef struct CropStruct
{
    /* Instance of the crop that is being planted */
    //int             cropUniqueIdentifier;
    char            cropName[128];

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
    
    int             userSeedingDate;
    int             userFloweringDate;
    int             userMaturityDate;
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

    double          calculatedFloweringTT;
    double          calculatedMaturityTT;
    double          calculatedSimAvgYield;
    double          calculatedSimMaxYield;
    double          calculatedSimMinYield;
    double          LWP_StressOnset;
    double          LWP_WiltingPoint;
    double          transpirationMax;

    int             harvestDateFinal;
    int             harvestCount;
    enum stage      stageGrowth;

    //int             rcCropNumber;
    //char            rcName[128];
    //int             rcActiveStatus;
    //int             rcYear;
    //int             rcDoy;
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
    //double          totalRealizedCrops;
} CropStruct;

typedef struct CommunityStruct
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
    //double          svUnstressedShootDailyGrowth;
    //double          svUnstressedRootDailyGrowth;
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

    CropStruct     *Crop;
    int             NumCrop;
    int             NumActiveCrop;
} CommunityStruct;

typedef struct FieldOperationStruct
{
    int             opYear;
    int             opDay;
    int             status;

    /* Planting Order */
    /* cropName and plantID are shared with forced harvest structure */
    char            cropName[128];
    int             usesAutoIrrigation;
    int             usesAutoFertilization;
    int             plantID;
    double          plantingDensity;

    /* Tillage */
    char            opToolName[MAXSTRING];
    double          opDepth;
    double          opSDR;
    double          opMixingEfficiency;
    char            cropNameT[128];
    double          fractionThermalTime;
    double          killEfficiency;
    int             grainHarvest;
    double          forageHarvest;

    /* Fixed Irrigation */
    double          opVolume;

    /* Fixed Fertilization */
    char            opSource[MAXSTRING];
    double          opMass;
    char            opForm[MAXSTRING];
    char            opMethod[MAXSTRING];
    int             opLayer;    /* Starting from 1 */
    double          opC_Organic;
    double          opC_Charcoal;
    double          opN_Organic;
    double          opN_Charcoal;
    double          opN_NH4;
    double          opN_NO3;
    double          opP_Organic;
    double          opP_Charcoal;
    double          opP_Inorganic;
    double          opK;
    double          opS;
} FieldOperationStruct;

typedef struct autoIrrigationStruct
{
    char            cropName[128];
    int             startDay;
    int             stopDay;
    double          waterDepletion;
    int             lastSoilLayer;
} autoIrrigationStruct;
//
//typedef struct autoFertilizationStruct
//{
//    char            cropName[128];
//    int             startDay;
//    int             stopDay;
//    double          mass;
//    char            source[MAXSTRING];
//    char            form[MAXSTRING];
//    char            method[MAXSTRING];
//} autoFertilizationStruct;

typedef struct CropManagementStruct
{
    FieldOperationStruct *FixedFertilization;
    int             numFertilization;

    FieldOperationStruct *FixedIrrigation;
    int             numIrrigation;

    FieldOperationStruct *Tillage;
    int             numTillage;
    double         *tillageFactor;

    FieldOperationStruct *plantingOrder;
    int             totalCropsPerRotation;

    FieldOperationStruct *ForcedHarvest;
    int             numHarvest;

    autoIrrigationStruct *autoIrrigation;
    int             usingAutoIrr;
    int             usingAutoFert;

    //char            nextCropName[128];
    //int             nextCropSeedingDate;
    //int             nextCropSeedingYear;
} CropManagementStruct;

typedef struct SnowStruct
{
    double          Snow;
    double          snowFall;
    double          snowMelt;
    double          snowCover;
    double          snowEvaporationVol;
} SnowStruct;

typedef struct ResidueStruct
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

    double         *residueAbgd;
    double         *residueRt;
    double         *residueRz;
    double         *residueAbgdN;
    double         *residueRtN;
    double         *residueRzN;
    double          yearResidueBiomass;
    double          yearResidueHarvested;
    double          yearRootBiomass;
    double          yearRhizodepositionBiomass;
    double         *manureC;
    double         *manureN;    /* Mg/ha */
} ResidueStruct;

typedef struct SoilCarbonStruct
{
    double         *factorComposite;
    double         *carbonRespired;
    double         *rootBiomassInput;
    double         *rhizBiomassInput;
    double         *abgdBiomassInput;
    double         *rootCarbonInput;
    double         *rhizCarbonInput;
    double         *manuCarbonInput;
    double         *abgdCarbonInput;
    double         *carbonMassInitial;
    double         *carbonMassFinal;
    double         *annualDecompositionFactor;
    double         *annualSoilCarbonDecompositionRate;
    double         *annualCarbonInputByLayer;
    double         *annualHumifiedCarbonMass;
    double         *annualRespiredCarbonMass;
    double         *annualRespiredResidueCarbonMass;
    double         *annualHumificationCoefficient;
    double         *annualNmineralization;
    double         *annualNImmobilization;
    double         *annualNNetMineralization;
    double          annualAmmoniumNitrification;
    double          annualNitrousOxidefromNitrification;
    double          annualAmmoniaVolatilization;
    double          annualNO3Denitrification;
    double          annualNitrousOxidefromDenitrification;
    double          annualNitrateLeaching;
    double          annualAmmoniumLeaching;
} SoilCarbonStruct;

typedef struct WeatherStruct
{
    double          siteAltitude;
    double          siteLatitude;
    double          screeningHeight;
    int             length;
    double         *yearlyAmplitude;
    double         *annualAverageTemperature;
    int            *lastDoy;
    double        **wind;
    double        **ETref;
    double        **precipitation;
    double        **RHmax;
    double        **RHmin;
    double        **solarRadiation;
    double        **tMax;
    double        **tMin;
    double          atmosphericPressure;
} WeatherStruct;

typedef struct PrintStruct
{
    char	    var_name[16];
    char	    unit[16];
    double	   *print_var;
} PrintStruct;

typedef struct SummaryStruct
{
    double          abgd_c_input;
    double          root_c_input;
    double          residue_biomass;
    double          produced_root;
    double          residue_resp;
    double          hum;
    double          soil_resp;
    double          n_mineralization;
    double          n_immobilization;
    double          n_net_mineralization;
    double          nh4_nitrification;
    double          n2o_from_nitrification;
    double          nh3_volatilization;
    double          no3_denirification;
    double          n2o_from_denitrification;
    double          no3_leaching;
    double          nh4_leaching;
    double          initial_soc;
    double          final_soc;
} SummaryStruct;

#ifdef _CYCLES_
typedef struct grid_struct
{
    SoilStruct      Soil;
    CropManagementStruct CropManagement;
    CommunityStruct Community;
    ResidueStruct   Residue;
    SoilCarbonStruct SoilCarbon;
    WeatherStruct   Weather;
    SnowStruct      Snow;
    SummaryStruct   Summary;
} grid_struct;
#endif

typedef struct CyclesStruct
{
    SimControlStruct SimControl;

#ifdef _CYCLES_
    grid_struct    *grid;
#else
    SoilStruct      Soil;
    CropManagementStruct CropManagement;
    CommunityStruct Community;
    ResidueStruct   Residue;
    SoilCarbonStruct SoilCarbon;
    WeatherStruct   Weather;
    SnowStruct      Snow;
    SummaryStruct   Summary;
#endif

    PrintStruct    *daily_output;
    PrintStruct    *annual_output;
} *CyclesStruct;

#endif
