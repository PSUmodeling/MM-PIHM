#ifndef ELEM_STRUCT_HEADER
#define ELEM_STRUCT_HEADER

/*****************************************************************************
 * Element attribute
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * soil_type                int         element soil type
 * lc_type                  int         element land cover type
 * bc_type                  int[]       element boundary condition type
 * meteo_type               int         element meteorological forcing type
 * lai_type                 int         element leaf area index forcing type
 ****************************************************************************/
typedef struct attrib_struct
{
    int             soil_type;
    int             lc_type;
    int             bc_type[NUM_EDGE];
    int             meteo_type;
    int             lai_type;
} attrib_struct;

/*****************************************************************************
 * Topography parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * area                     double      area of element [m2]
 * x                        double      x of centroid [m]
 * y                        double      y of centroid [m]
 * zmin                     double      bedrock elevation [m]
 * zmax                     double      surface elevation [m]
 * edge                     double[]    length of edge (Edge i is from node i
 *                                        to node i + 1) [m]
 * nabrdist_x               double[]    distance to neighbor in x direction
 *                                        [m]
 * nabrdist_y               double[]    distance to neighbor in y direction
 *                                        [m]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * slope                    double      slope of element [degree]
 * aspect                   double      surface aspect of element [degree]
 * svf                      double      sky view factor [-]
 * h_phi                    double[]    unobstructed angle in each direction
 *                                        [degree]
 * ---------------------------------------------------------------------------
 * Variables below only used in RT-Flux-PIHM
 * ---------------------------------------------------------------------------
 * areasub                  double[]
 ****************************************************************************/
typedef struct topo_struct
{
    double          area;
    double          x;
    double          y;
    double          zmin;
    double          zmax;
    double          edge[NUM_EDGE];
    double          nabrdist_x[NUM_EDGE];
    double          nabrdist_y[NUM_EDGE];
#ifdef _NOAH_
    double          slope;
    double          aspect;
    double          svf;
    double          h_phi[36];
#endif
#ifdef _RT_
    double          areasub[NUM_EDGE];
#endif
} topo_struct;

/*****************************************************************************
 * Soil parameter
 * ---------------------------------------------------------------------------
 *
 * Variables                Type        Description
 * ==========               ==========  ====================
 * depth                    double      soil depth [m]
 * ksath                    double      horizontal saturated hydraulic
 *                                        conductivity [m s-1]
 * ksatv                    double      vertical saturated hydraulic
 *                                        conductivity [m s-1]
 * kinfv                    double      saturated infiltration conductivity
 *                                        [m s-1]
 * dinf                     double      depth from ground surface accross which
 *                                        head gradient is calculated for
 *                                        infiltration [m]
 * alpha                    double      alpha from van Genuchten eqn [m-1]
 * beta                     double      beta (n) from van Genuchten eqn [-]
 * porosity                 double      soil porosity [m3 m-3]
 * smcmax                   double      maximum soil moisture content [m3 m-3]
 * smcmin                   double      residual soil moisture content
 *                                        [m3 m-3]
 * smcwlt                   double      wilting point [m3 m-3]
 * smcref                   double      soil moisture threshold where
 *                                        transpiration begins to stress
 *                                        [m3 m-3]
 * dmac                     double      macropore depth [m]
 * kmach                    double      macropore horizontal saturated
 *                                        hydraulic conductivity [m s-1]
 * kmacv                    double      macropore vertical saturated hydraulic
 *                                        conductivity [m s-1]
 * areafv                   double      macropore area fraction on a vertical
 *                                        cross-section [m2 m-2]
 * areafh                   double      macropore area fraction on a
 *                                        horizontal cross-section [m2 m-2]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * csoil                    double      soil heat capacity [J m-3 K-1]
 * quartz                   double      soil quartz content [-]
 * smcdry                   double      dry soil moisture threshold where
 *                                        direct evap frm top layer ends
 *                                        [m3 m-3]
 ****************************************************************************/
typedef struct soil_struct
{
    double          depth;
    double          ksath;
    double          ksatv;
    double          kinfv;
    double          dinf;
    double          alpha;
    double          beta;
    double          porosity;
    double          smcmax;
    double          smcmin;
    double          smcwlt;
    double          smcref;
    double          dmac;
    double          kmach;
    double          kmacv;
    double          areafv;
    double          areafh;
#ifdef _NOAH_
    double          csoil;
    double          quartz;
    double          smcdry;
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

/*****************************************************************************
 * Land cover parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * shdfac                   double      areal fractional coverage of green
 *                                        vegetation (0.0-1.0) [-]
 * shdmin                   double      minimum areal fractional coverage of
 *                                        green vegetation [-]
 * shdmax                   double      maximum areal fractional coverage of
 *                                        green vegetation [-]
 * laimin                   double      minimum LAI across all seasons for a
 *                                        vegetation type [m2 m-2]
 * laimax                   double      maximum LAI across all seasons for a
 *                                        vegetation type [m2 m-2]
 * snup                     double      threshold snow depth (in water
 *                                        equivalent) that implies 100% snow
 *                                        cover [m]
 * cfactr                   double      parameter used in the canopy
 *                                        inteception calculation [-]
 * emissmin                 double      minimum emissivity [-]
 * emissmax                 double      maximum emissivity [-]
 * albedomin                double      minimum background albedo [-]
 * albedomax                double      maximum background albedo [-]
 * z0min                    double      minimum roughness length [m]
 * z0max                    double      maximum roughness length [m]
 * rought                   double      surface roughness (Manning's n)
 *                                        [s m-1/3]
 * cmcfactr                 double      canopy water capacity per LAI [m]
 * bare                     int         flag that indicates bare ground
 * isurban                  int         flag that indicates urban
 ****************************************************************************/
typedef struct lc_struct
{
    double          shdfac;
    double          shdmin;
    double          shdmax;
    double          laimin;
    double          laimax;
    double          snup;
    double          cfactr;
    double          emissmax;
    double          emissmin;
    double          albedomax;
    double          albedomin;
    double          z0max;
    double          z0min;
    double          rough;
    double          cmcfactr;
    int             bare;
    int             isurban;
} lc_struct;

/*****************************************************************************
 * Ecophysiological parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * rsmin                    double      minimum canopy resistance [s m-1]
 * rgl                      double      reference incoming solar flux for
 *                                        photosynthetically active canopy
 *                                        [W m-2]
 * hs                       double      parameter used in vapor pressure
 *                                        deficit function [-]
 * topt                     double      optimum transpiration air temperature
 *                                        [K]
 * rsmax                    double      cuticular resistance [s m-1]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM-BGC
 * ---------------------------------------------------------------------------
 * woody                    int         flag: 1 = woody, 0 = non-woody
 * evergreen                int         flag: 1 = evergreen, 0 = deciduous
 * c3_flag                  int         flag: 1 = C3,  0 = C4
 * phenology_flag           int         flag: 1 = phenology model, 0 = user
 *                                        defined
 * onday                    int         day of year when leaves on
 * offday                   int         day of year yearday leaves off
 * transfer_days            int         growth period for transfer [day]
 * litfall_days             int         growth period for litfall [day]
 * leaf_turnover            double      annual leaf turnover fraction [yr-1]
 * froot_turnover           double      annual fine root turnover fraction
 *                                        [yr-1]
 * livewood_turnover        double      annual live wood turnover fraction
 *                                        [yr-1]
 * daily_mortality_turnover double      daily mortality turnover [day-1]
 * daily_fire_turnover      double      daily fire turnover [day-1]
 * alloc_frootc_leafc       double      new fine root C to new leaf C [-]
 * alloc_newstemc_newleafc  double      new stem C to new leaf C [-]
 * alloc_newlivewoodc_newwoodc
 *                          double      new livewood C:new wood C [-]
 * alloc_crootc_stemc       double      new live croot C to new live stem C
 *                                        [-]
 * alloc_prop_curgrowth     double      daily allocation to current growth [-]
 * avg_proj_sla             double      canopy average projected SLA
 *                                        [m2 kgC-1]
 * sla_ratio                double      ratio of shaded to sunlit projected
 *                                        SLA [-]
 * lai_ratio                double      ratio of (all-sided LA / one-sided LA)
 *                                        [-]
 * ext_coef                 double      canopy light extinction coefficient
 *                                        [-]
 * flnr                     double      leaf N in Rubisco [kgNRub kgNleaf-1]
 * psi_open                 double      psi at start of conductance reduction
 *                                        [MPa]
 * psi_close                double      psi at complete conductance reduction
 *                                        [MPa]
 * vpd_open                 double      vpd at start of conductance reduction
 *                                        [Pa]
 * vpd_close                double      vpd at complete conductance reduction
 *                                        [Pa]
 * froot_cn                 double      C:N for fine roots [kgC kgN-1]
 * leaf_cn                  double      C:N for leaves [kgC kgN-1]
 * livewood_cn              double      C:N for live wood [kgC kgN-1]
 * deadwood_cn              double      C:N for dead wood [kgC kgN-1]
 * leaflitr_cn              double      constant C:N for leaf litter
 *                                        [kgC kgN-1]
 * leaflitr_flab            double      leaf litter labile fraction [-]
 * leaflitr_fucel           double      leaf litter unshielded cellulose
 *                                        fraction [-]
 * leaflitr_fscel           double      leaf litter shielded cellulose
 *                                        fraction [-]
 * leaflitr_flig            double      leaf litter lignin fraction [-]
 * frootlitr_flab           double      fine root litter labile fraction [-]
 * frootlitr_fucel          double      fine root litter unshielded cellulose
 *                                        fraction [-]
 * frootlitr_fscel          double      fine root litter shielded cellulose
 *                                        fraction [-]
 * frootlitr_flig           double      fine root litter lignin fraction [-]
 * deadwood_fucel           double      dead wood unshileded cellulose
 *                                        fraction [-]
 * deadwood_fscel           double      dead wood shielded cellulose fraction
 *                                        [-]
 * deadwood_flig            double      dead wood lignin fraction [-]
 ****************************************************************************/
typedef struct epconst_struct
{
    double          rsmin;
    double          rgl;
    double          hs;
    double          topt;
    double          rsmax;
#ifdef _BGC_
    int             woody;
    int             evergreen;
    int             c3_flag;
    int             phenology_flag;
    int             onday;
    int             offday;
    int             transfer_days;
    int             litfall_days;
    double          leaf_turnover;
    double          froot_turnover;
    double          livewood_turnover;
    double          daily_mortality_turnover;
    double          daily_fire_turnover;
    double          alloc_frootc_leafc;
    double          alloc_newstemc_newleafc;
    double          alloc_newlivewoodc_newwoodc;
    double          alloc_crootc_stemc;
    double          alloc_prop_curgrowth;
    double          avg_proj_sla;
    double          sla_ratio;
    double          lai_ratio;
    double          ext_coef;
    double          flnr;
    double          psi_open;
    double          psi_close;
    double          vpd_open;
    double          vpd_close;
    double          froot_cn;
    double          leaf_cn;
    double          livewood_cn;
    double          deadwood_cn;
    double          leaflitr_cn;
    double          leaflitr_flab;
    double          leaflitr_fucel;
    double          leaflitr_fscel;
    double          leaflitr_flig;
    double          frootlitr_flab;
    double          frootlitr_fucel;
    double          frootlitr_fscel;
    double          frootlitr_flig;
    double          deadwood_fucel;
    double          deadwood_fscel;
    double          deadwood_flig;
#endif
} epconst_struct;

/*****************************************************************************
 * Physical states
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * rzd                      double      rooting depth [m]
 * macpore_status           int         macropore status
 * rc                       double      canopy resistance [s m-1]
 * pc                       double      plant coefficient [-]
 * proj_lai                 double      live projected leaf area index
 *                                        [m2 m-2]
 * rcs                      double      incoming solar rc factor [-]
 * rct                      double      air temperature rc factor [-]
 * rcq                      double      vapor pressure deficit rc factor [-]
 * rcsoil                   double      soil moisture rc factor [-]
 * albedo                   double      surface albedo including snow effect
 *                                        [-]
 * zlvl                     double      height above ground of atmospheric
 *                                        forcing variables [m]
 * zlvl_wind                double      height above ground of wind
 *                                        observations [m]
 * sfcspd                   double      wind speed at height zlvl above ground
 *                                        [m s-1]
 * rh                       double      relative humidity [100%]
 * sfcprs                   double      surface pressure at height zlvl abouve
 *                                        ground [Pa]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * alb                      double      backround snow-free surface albedo [-]
 * snoalb                   double      upper bound on maximum albedo over
 *                                        deep snow [-]
 * nroot                    int         number of root layers, a function of
 *                                        vegetation type
 * rtdis                    double[]    root distribution [-]
 * nsoil                    int         number of soil layers
 * sldpth                   double[]    thickness of each soil layer [m]
 * soilw                    double      available soil moisture in root zone
 *                                        (fraction between smcwlt and smcmax)
 *                                        [-]
 * frzk                     double      frozen ground parameter [-]
 * frzx                     double      adjusted frozen ground parameter [-]
 * czil                     double      zilitinkevich constant [-]
 * emissi                   double      surface emissivity (between 0 and 1)
 *                                        [-]
 * ch                       double      surface exchange coefficient for heat
 *                                        and moisture (m s-1)
 * cm                       double      surface exchange coefficient for
 *                                        momentum [m s-1]
 * rch                      double      = ch * air density * CP [W m-2 K-1]
 * z0                       double      time varying roughness length as
 *                                        function of snow depth [-]
 * fcr                      double      reduction of infiltration caused by
 *                                        frozen ground [-]
 * nmacd                    int         number of soil layers with macropore
 * salp                     double      shape parameter of distribution
 *                                        function of snow cover [-]
 * fxexp                    double      soil evaporation exponent used in
 *                                        direct evaporation [-]
 * sbeta                    double      parameter used to calculate vegetation
 *                                        effect on soil heat [-]
 * lvcoef                   double      parameter controls surface snow albedo
 *                                        in the presence of snowcover [-]
 * snotime1                 double      age of the snow on the ground [s]
 * ribb                     double      bulk richardson number used to limit
 *                                        the dew/frost [-]
 * beta                     double      ratio of actual/potential evap [-]
 * sncovr                   double      fractional snow cover [-]
 * q1                       double      effective mixing ratio at surface
 *                                        [kg kg-1]
 * q2                       double      mixing ratio at height zlvl above
 *                                        [kg kg-1]
 * ffrozp                   double      fraction of frozen precipitation [-]
 * z0brd                    double      background fixed roughness length [-]
 * embrd                    double      background surface emissivity [-]
 * q2sat                    double      saturation air humidity at height zlvl
 *                                        above ground [kg kg-1]
 * q2d                      double      air humidity deficit [kg kg-1]
 * dqsdt2                   double      slope of saturation specific humidity
 *                                        curve at T = sfctmp [kg kg-1 K-1]
 * nwtbl                    int         layer where water table is within
 * sndens                   double      snow density (dimensionless fraction
 *                                        of H2O density) [-]
 * snowh                    double      actual snow depth [-]
 * sncond                   double      snow thermal conductivity [W m-1 K-1]
 * rr                       double      parameter in Penman potential
 *                                        evaopration [-]
 * epsca                    double      parameter in Penman potential
 *                                        evaopration [K]
 * eta_kinematic            double      atctual latent heat flux [kg m-2 s-1]
 * zbot                     double      depth of lower boundary soil
 *                                        temperature [m]
 * tbot                     double      bottom soil temperature (local yearly-
 *                                        mean sfc air temperature) [K]
 * gwet                     double      fraction of transpiration from
 *                                        groundwater [-]
 * satdpth                  double[]    depth of groundwater in each soil
 *                                        layer [m]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM-BGC
 * ---------------------------------------------------------------------------
 * co2                      double      atmospheric CO2 concentration [ppm]
 * ppfd_per_plaisun         double      ppfd per unit sunlit proj LAI
 *                                        [umol m-2 s-1]
 * ppfd_per_plaishade       double      ppfd per unit shaded proj LAI
 *                                        [umol m-2 s-1]
 * all_lai                  double      live all-sided leaf area index
 *                                        [m2 m-2]
 * plaisun                  double      sunlit projected leaf area index
 *                                        [m2 m-2]
 * plaishade                double      shaded projected leaf area index
 *                                        [m2 m-2]
 ****************************************************************************/
typedef struct pstate_struct
{
    double          rzd;
    int             macpore_status;
    double          rc;
    double          pc;
    double          proj_lai;
    double          rcs;
    double          rct;
    double          rcq;
    double          rcsoil;
    double          albedo;
    double          zlvl;
    double          zlvl_wind;
    double          sfcspd;
    double          rh;
    double          sfcprs;
#ifdef _NOAH_
    double          alb;
    double          snoalb;
    int             nroot;
    double          rtdis[MAXLYR];
    int             nsoil;
    double          sldpth[MAXLYR];
    double          soilw;
    double          frzk;
    double          frzx;
    double          czil;
    double          emissi;
    double          ch;
    double          rch;
    double          cm;
    double          z0;
    double          fcr;
    int             nmacd;
    double          salp;
    double          fxexp;
    double          sbeta;
    double          lvcoef;
    double          snotime1;
    double          ribb;
    double          beta;
    double          sncovr;
    double          q1;
    double          q2;
    double          ffrozp;
    double          z0brd;
    double          embrd;
    double          q2sat;
    double          q2d;
    double          dqsdt2;
    int             nwtbl;
    double          sndens;
    double          snowh;
    double          sncond;
    double          rr;
    double          epsca;
    double          eta_kinematic;
    double          zbot;
    double          tbot;
    double          gwet;
    double          satdpth[MAXLYR];
#endif
#ifdef _BGC_
    double          co2;
    double          ppfd_per_plaisun;
    double          ppfd_per_plaishade;
    double          all_lai;
    double          plaisun;
    double          plaishade;
#endif
} pstate_struct;

/*****************************************************************************
 * Water states
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * surf                     double      surface water level [m]
 * unsat                    double      unsaturated zone water storage [m]
 * gw                       double      groundwater level [m]
 * sneqv                    double      liquid water-equivalent snow depth [m]
 * cmcmax                   double      maximum canopy water capacity [m]
 * cmc                      double      interception storage [m]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * smc                      double[]    total soil moisture content [m3 m-3]
 * sh2o                     double[]    unfrozen soil moisture content
 *                                        [m3 m-3]
 * soilm                    double      total soil column moisture content [m]
 ****************************************************************************/
typedef struct wstate_struct
{
    double          surf;
    double          unsat;
    double          gw;
    double          sneqv;
    double          cmcmax;
    double          cmc;
#ifdef _NOAH_
    double          smc[MAXLYR];
    double          sh2o[MAXLYR];
    double          soilm;
#endif
} wstate_struct;

/*****************************************************************************
 * Water fluxes
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * ovlflow                  double[]    overland flow [m3 s-1]
 * subsurf                  double[]    subsurface flow [m3 s-1]
 * prcp                     double      precepitation on each element [m s-1]
 * pcpdrp                   double      combined prcp and drip (from canopy)
 *                                        that goes into the soil [m s-1]
 * infil                    double      variable infiltration rate [m s-1]
 * rechg                    double      recharge rate to groundwater [m s-1]
 * drip                     double      through-fall of precipitation and/or
 *                                        dew [m s-1]
 * edir                     double      direct soil evaporation [m s-1]
 * ett                      double      total plant transpiration [m s-1]
 * ec                       double      canopy water evaporation [m s-1]
 * etp                      double      potential evaporation [m s-1]
 * eta                      double      actual evapotranspiration [m s-1]
 * edir_surf                double      direct evaporation from surface water
 *                                        [m s-1]
 * edir_unsat               double      direct evaporation from unsaturated
 *                                        zone [m s-1]
 * edir_gw                  double      direct evaporation from saturated zone
 *                                        [m s-1]
 * ett_unsat                double      transpiration from unsaturated zone
 *                                        [m s-1]
 * ett_gw                   double      transpiration from saturated zone
 *                                        [m s-1]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * et                       doube[]     plant transpiration from each soil
 *                                        layer [m s-1]
 * runoff2                  double      total subsurface flow [m s-1]
 * runoff2_lyr              double[]    subsurface flow from each soil layer
 *                                        [m s-1]
 * runoff3                  double      numerical trunctation in excess of
 *                                        porosity (smcmax) for a given soil
 *                                        layer at the end of a time step
 *                                        [m s-1]
 * smflxv                   double[]    vertical soil mositure flux between
 *                                        soil layers [m s-1]
 * smflxh                   double[][]  horizontal soil moisture flux at each
 *                                        soil layer from each edge [m s-1]
 * dew                      double      dewfall (or frostfall for T < 273.15)
 *                                        [m s-1]
 * snomlt                   double      water equivalent snow melt [m s-1]
 * esnow                    double      sublimation from (or deposition to)
 *                                        snowpack [m s-1];
 * etns                     double      [m s-1]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM-Cycles
 * ---------------------------------------------------------------------------
 * eres                     double      evaporation from residue [m s-1]
 ****************************************************************************/
typedef struct wflux_struct
{
    double          ovlflow[NUM_EDGE];
    double          subsurf[NUM_EDGE];
    double          prcp;
    double          pcpdrp;
    double          infil;
    double          rechg;
    double          drip;
    double          edir;
    double          ett;
    double          ec;
    double          etp;
    double          eta;
    double          edir_surf;
    double          edir_unsat;
    double          edir_gw;
    double          ett_unsat;
    double          ett_gw;
#ifdef _NOAH_
    double          et[MAXLYR];
    double          runoff2;
    double          runoff2_lyr[MAXLYR];
    double          runoff3;
    double          smflxv[MAXLYR];
    double          smflxh[NUM_EDGE][MAXLYR];
    double          dew;
    double          snomlt;
    double          esnow;
    double          etns;
#endif
#ifdef _CYCLES_
    double          eres;
#endif
} wflux_struct;

/*****************************************************************************
 * Water states
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * sfctmp                   double      air temperature at height zlvl above
 *                                        ground [K]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * t1                       double      ground/canopy/snowpack effective skin
 *                                        temperature [K]
 * th2                      double      air potential temperature at height
 *                                        zlvl above ground [K]
 * stc                      double[]    soil temperature [K]
 ****************************************************************************/
typedef struct estate_struct
{
    double          sfctmp;
#ifdef _NOAH_
    double          t1;
    double          th2;
    double          stc[MAXLYR];
#endif
} estate_struct;

/*****************************************************************************
 * Water states
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * soldn                    double      solar downward radiation [W m-2]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * solnet                   double      net downward solar radiation [W m-2]
 * etp                      double      potential evaporation [W m-2]
 * ssoil                    double      soil heat flux [W m-2]
 * eta                      double      actual latent heat flux [W m-2]
 * sheat                    double      sensible heat flux [W m-2]
 * fdown                    double      radiation forcing at the surface
 *                                        [W m-2]
 * lwdn                     double      absorbed longwave downward radiation
 *                                        [W m-2]
 * ec                       double      canopy water evaporation [W m-2]
 * edir                     double      direct soil evaporation [W m-2]
 * et                       double[]    plant transpiration from each soil
 *                                        layer [W m-2]
 * ett                      double      total plant transpiration [W m-2]
 * esnow                    double      sublimation from (or deposition)
 *                                        snowpack [W m-2]
 * soldir                   double      direct solar radiation [W m-2]
 * soldif                   double      diffused solar radiation [W m-2]
 * longwave                 double      longwave radiation forcing [W m-2]
 * flx1                     double      latent heat flux from precipitation
 *                                        accumulating as snow [W m-2]
 * flx2                     double      freezing rain latent heat flux [W m-2]
 * flx3                     double      snow melt latent heat flux [W m-2]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM-BGC
 * ---------------------------------------------------------------------------
 * swabs                    double      canopy absorbed shortwave flux [W m-2]
 * swtrans                  double      transmitted shortwave flux [W m-2]
 * swabs_per_plaisun        double      swabs per unit sunlit proj LAI [W m-2]
 * swabs_per_plaishade      double      swabs per unit shaded proj LAI [W m-2]
 * parabs                   double      PAR absorbed by canopy [W m-2]
 ****************************************************************************/
typedef struct eflux_struct
{
    double          soldn;
#ifdef _NOAH_
    double          solnet;
    double          etp;
    double          ssoil;
    double          eta;
    double          sheat;
    double          fdown;
    double          lwdn;
    double          ec;
    double          edir;
    double          et[MAXLYR];
    double          ett;
    double          esnow;
    double          soldir;
    double          soldif;
    double          longwave;
    double          flx1;
    double          flx2;
    double          flx3;
#endif
#ifdef _BGC_
    double          swabs;
    double          swtrans;
    double          swabs_per_plaisun;
    double          swabs_per_plaishade;
    double          parabs;
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
    double          svRootingDepth;     /* maximum */
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

/*****************************************************************************
 * Boundary conditions
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * head                     double[]    value of Dirichlet-type boudnary
 *                                        condition [m]
 * flux                     double[]    value of Neumann-type boundary
 *                                        condition [m3 s-1]
 ****************************************************************************/
typedef struct bc_struct
{
    double          head[NUM_EDGE];
    double          flux[NUM_EDGE];
} bc_struct;

/*****************************************************************************
 * Initial conditions
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * cmc                      double      interception storage [m]
 * sneqv                    double      liquid water-equivalent snow depth [m]
 * surf                     double      surface water level [m]
 * unsat                    double      unsaturated zone water storage [m]
 * gw                       double      groundwater level [m]
 * ---------------------------------------------------------------------------
 * Variables below only used in Flux-PIHM
 * ---------------------------------------------------------------------------
 * t1                       double      ground/canopy/snowpack effective skin
 *                                        temperature [K]
 * snowh                    double      actual snow depth [-]
 * stc                      double[]    soil temperature [K]
 * smc                      double[]    total soil moisture content [m3 m-3]
 * sh2o                     double[]    unfrozen soil moisture content
 *                                        [m3 m-3]
 ****************************************************************************/
typedef struct ic_struct
{
    double          cmc;
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
} ic_struct;

#ifdef _BGC_
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
} bgcic_struct;
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
    double          soluteFluxLat[4][MAXLYR];
    double          soluteFluxVert[MAXLYR];
} solute_struct;
#endif

#ifdef _DAILY_
/*****************************************************************************
 * Daily average
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * counter                  int         counter used for averaging
 * daylight_counter         int         counter used for daytime averaging
 * avg_surf                 double      daily average surface water level [m]
 * avg_unsat                double      daily average unsaturated water
 *                                        storage [m]
 * avg_gw                   double      daily average groundwater level [m]
 * avg_sh2o                 double      daily average unfrozen soil moisture
 *                                        content [m3 m-3]
 * avg_ovlflow              double[]    daily average overland flow [m3 s-1]
 * avg_subsurf              double[]    daily average subsurface flow [m3 s-1]
 * avg_et                   double[]    daily average evaportranspiration
 *                                        [m s-1]
 * avg_smflxv               double[]    daily average vertical soil moisture
 *                                        flux [m s-1]
 * dayl                     double      day length [s]
 * prev_dayl                double      day length of previous day [s]
 * avg_q2d                  double      daily average mixing ratio deficit
 *                                        [kg kg-1]
 * avg_sfcprs               double      daily average air pressure [Pa]
 * avg_ch                   double      daily average surface exchange
 *                                        coefficient for heat and moisture
 *                                        [m s-1]
 * avg_albedo               double      daily average surface albedo [-]
 * avg_sfcspd               double      daily average wind speed [m s-1]
 * avg_sncovr               double      daily average snow cover fraction [-]
 * tmax                     double      daily maximum air temperature [K]
 * tmin                     double      daily minimum air temperature [K]
 * avg_sfctmp               double      daily average air temperature [K]
 * tday                     double      daytime average air temperature [K]
 * tnight                   double      nighttime average air temperature [K]
 * avg_stc[MAXLYR]          double      daily average soil temperature [K]
 * avg_soldn                double      daytime average downward solar
 *                                        radiation [W m-2]
 * solar_total              double      daily total solar radiation [J m-2]
 ****************************************************************************/
typedef struct daily_struct
{
    int             counter;
    int             daylight_counter;

    double          avg_surf;
    double          avg_unsat;
    double          avg_gw;
    double          avg_sh2o[MAXLYR];

    double          avg_ovlflow[NUM_EDGE];
    double          avg_subsurf[NUM_EDGE];
    double          avg_et[MAXLYR];
    double          avg_smflxv[MAXLYR];

    double          dayl;
    double          prev_dayl;
    double          avg_q2d;
    double          avg_sfcprs;
    double          avg_ch;
    double          avg_albedo;
    double          avg_sfcspd;
    double          avg_sncovr;

    double          tmax;
    double          tmin;
    double          avg_sfctmp;
    double          tday;
    double          tnight;
    double          avg_stc[MAXLYR];

    double          avg_soldn;
    double          solar_total;
} daily_struct;
#endif

#ifdef _BGC_
/*****************************************************************************
 * Storage of daily average (only used in Flux-PIHM-BGC)
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * flag                     int*        flag indicating storage status
 * dayl                     double*     day length [s]
 * prev_dayl                double*     day length of previous day [s]
 * tmax                     double*     daily maximum air temperature [K]
 * tmin                     double*     daily minimum air temperature [K]
 * sfctmp                   double*     daily average air temperature [K]
 * tday                     double*     daytime average air temperature [K]
 * tnight                   double*     nighttime average air temperature [K]
 * q2d                      double*     daily average mixing ratio deficit
 *                                *       [kg kg-1]
 * sfcprs                   double*     daily average air pressure [Pa]
 * avg_soldn                double*     daytime average downward solar
 *                                *       radiation [W m-2]
 * avg_stc[MAXLYR]          double*     daily average soil temperature [K]
 * avg_sh2o                 double*     daily average unfrozen soil moisture
 *                                *       content [m3 m-3]
 * avg_surf                 double*     daily average surface water level [m]
 * avg_unsat                double*     daily average unsaturated water
 *                                *       storage [m]
 * avg_gw                   double*     daily average groundwater level [m]
 * avg_albedo               double*     daily average surface albedo [-]
 * avg_ch                   double*     daily average surface exchange
 *                                       coefficient for heat and moisture
 *                                       [m s-1]
 * surfflx                  double[]*   daily average overland flow [m3 s-1]
 * subsurfflx               double[]*   daily average subsurface flow [m3 s-1]
 ****************************************************************************/
typedef struct stor_struct
{
    int            *flag;
    double         *dayl;       /* (s)     daylength */
    double         *prev_dayl;
    double         *tmax;       /* (deg C) daily maximum air temperature */
    double         *tmin;       /* (deg C) daily minimum air temperature */
    double         *sfctmp;     /* (deg C) daily average temperature */
    double         *tday;
    double         *tnight;
    double         *q2d;        /* (m3/m3) mixing ratio deficit */
    double         *sfcprs;
    double         *soldn;      /* (W/m2)  daylight avg shortwave flux density */
    double         *stc[MAXLYR];
    double         *sh2o[MAXLYR];
    double         *surf;
    double         *unsat;
    double         *gw;
    double         *albedo;
    double         *ch;
    double         *surfflx[NUM_EDGE];
    double         *subsurfflx[NUM_EDGE];
} stor_struct;


/*****************************************************************************
 * Carbon state variables (including sums for sources and sinks)
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * leafc                    double      leaf C [kgC m-2]
 * leafc_storage            double      leaf C storage [kgC m-2]
 * leafc_transfer           double      leaf C transfer [kgC m-2]
 * frootc                   double      fine root C [kgC m-2]
 * frootc_storage           double      fine root C storage [kgC m-2]
 * frootc_transfer          double      fine root C transfer [kgC m-2]
 * livestemc                double      live stem C [kgC m-2]
 * livestemc_storage        double      live stem C storage [kgC m-2]
 * livestemc_transfer       double      live stem C transfer [kgC m-2]
 * deadstemc                double      dead stem C [kgC m-2]
 * deadstemc_storage        double      dead stem C storage [kgC m-2]
 * deadstemc_transfer       double      dead stem C transfer [kgC m-2]
 * livecrootc               double      live coarse root C [kgC m-2]
 * livecrootc_storage       double      live coarse root C storage [kgC m-2]
 * livecrootc_transfer      double      live coarse root C transfer [kgC m-2]
 * deadcrootc               double      dead coarse root C [kgC m-2]
 * deadcrootc_storage       double      dead coarse root C storage [kgC m-2]
 * deadcrootc_transfer      double      dead coarse root C transfer [kgC m-2]
 * gresp_storage            double      growth respiration storage [kgC m-2]
 * gresp_transfer           double      growth respiration transfer [kgC m-2]
 * cwdc                     double      coarse woody debris C [kgC m-2]
 * litr1c                   double      litter labile C [kgC m-2]
 * litr2c                   double      litter unshielded cellulose C [kgC m-2]
 * litr3c                   double      litter shielded cellulose C [kgC m-2]
 * litr4c                   double      litter lignin C [kgC m-2]
 * soil1c                   double      microbial recycling pool C (fast) [kgC m-2]
 * soil2c                   double      microbial recycling pool C (medium) [kgC m-2]
 * soil3c                   double      microbial recycling pool C (slow) [kgC m-2]
 * soil4c                   double      recalcitrant SOM C (humus, slowest) [kgC m-2]
 * cpool                    double      temporary photosynthate C pool [kgC m-2]
 * psnsun_src               double      SUM of gross PSN from sulit canopy [kgC m-2]
 * psnshade_src             double      SUM of gross PSN from shaded canopy [kgC m-2]
 * leaf_mr_snk              double      SUM of leaf maint resp [kgC m-2]
 * leaf_gr_snk              double      SUM of leaf growth resp [kgC m-2]
 * froot_mr_snk             double      SUM of fine root maint resp [kgC m-2]
 * froot_gr_snk             double      SUM of fine root growth resp [kgC m-2]
 * livestem_mr_snk          double      SUM of live stem maint resp [kgC m-2]
 * livestem_gr_snk          double      SUM of live stem growth resp [kgC m-2]
 * deadstem_gr_snk          double      SUM of dead stem growth resp [kgC m-2]
 * livecroot_mr_snk         double      SUM of live coarse root maint resp [kgC m-2]
 * livecroot_gr_snk         double      SUM of live coarse root growth resp [kgC m-2]
 * deadcroot_gr_snk         double      SUM of dead coarse root growth resp [kgC m-2]
 * litr1_hr_snk             double      SUM of labile litr microbial resp [kgC m-2]
 * litr2_hr_snk             double      SUM of cellulose litr microbial resp [kgC m-2]
 * litr4_hr_snk             double      SUM of lignin litr microbial resp [kgC m-2]
 * soil1_hr_snk             double      SUM of fast microbial respiration [kgC m-2]
 * soil2_hr_snk             double      SUM of medium microbial respiration [kgC m-2]
 * soil3_hr_snk             double      SUM of slow microbial respiration [kgC m-2]
 * soil4_hr_snk             double      SUM of recalcitrant SOM respiration [kgC m-2]
 * fire_snk                 double      SUM of fire losses [kgC m-2]
 ****************************************************************************/
typedef struct cstate_struct
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
    double          psnsun_src;
    double          psnshade_src;
    double          leaf_mr_snk;
    double          leaf_gr_snk;
    double          froot_mr_snk;
    double          froot_gr_snk;
    double          livestem_mr_snk;
    double          livestem_gr_snk;
    double          deadstem_gr_snk;
    double          livecroot_mr_snk;
    double          livecroot_gr_snk;
    double          deadcroot_gr_snk;
    double          litr1_hr_snk;
    double          litr2_hr_snk;
    double          litr4_hr_snk;
    double          soil1_hr_snk;
    double          soil2_hr_snk;
    double          soil3_hr_snk;
    double          soil4_hr_snk;
    double          fire_snk;
} cstate_struct;

/*****************************************************************************
 * Daily carbon flux variables
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * m_leafc_to_litr1c        double      mortality flux [kgC m-2 day-1]
 * m_leafc_to_litr2c        double      mortality flux [kgC m-2 day-1]
 * m_leafc_to_litr3c        double      mortality flux [kgC m-2 day-1]
 * m_leafc_to_litr4c        double      mortality flux [kgC m-2 day-1]
 * m_frootc_to_litr1c       double      mortality flux [kgC m-2 day-1]
 * m_frootc_to_litr2c       double      mortality flux [kgC m-2 day-1]
 * m_frootc_to_litr3c       double      mortality flux [kgC m-2 day-1]
 * m_frootc_to_litr4c       double      mortality flux [kgC m-2 day-1]
 * m_leafc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_frootc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livestemc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadstemc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livecrootc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadcrootc_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_leafc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_frootc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livestemc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadstemc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livecrootc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadcrootc_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livestemc_to_cwdc
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadstemc_to_cwdc
 *                          double      mortality flux [kgC m-2 day-1]
 * m_livecrootc_to_cwdc
 *                          double      mortality flux [kgC m-2 day-1]
 * m_deadcrootc_to_cwdc     double      mortality flux [kgC m-2 day-1]
 * m_gresp_storage_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_gresp_transfer_to_litr1c
 *                          double      mortality flux [kgC m-2 day-1]
 * m_leafc_to_fire          double      fire flux [kgC m-2 day-1]
 * m_frootc_to_fire         double      fire flux [kgC m-2 day-1]
 * m_leafc_storage_to_fire  double      fire flux [kgC m-2 day-1]
 * m_frootc_storage_to_fire double      fire flux [kgC m-2 day-1]
 * m_livestemc_storage_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_deadstemc_storage_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_livecrootc_storage_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_deadcrootc_storage_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_leafc_transfer_to_fire double      fire flux [kgC m-2 day-1]
 * m_frootc_transfer_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_livestemc_transfer_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_deadstemc_transfer_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_livecrootc_transfer_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_deadcrootc_transfer_to_fire
 *                          double      fire flux [kgC m-2 day-1]
 * m_livestemc_to_fire      double      fire flux [kgC m-2 day-1]
 * m_deadstemc_to_fire      double      fire flux [kgC m-2 day-1]
 * m_livecrootc_to_fire     double      fire flux [kgC m-2 day-1]
 * m_deadcrootc_to_fire     double      fire flux [kgC m-2 day-1]
 * m_gresp_storage_to_fire  double      fire flux [kgC m-2 day-1]
 * m_gresp_transfer_to_fire double      fire flux [kgC m-2 day-1]
 * m_litr1c_to_fire         double      fire flux [kgC m-2 day-1]
 * m_litr2c_to_fire         double      fire flux [kgC m-2 day-1]
 * m_litr3c_to_fire         double      fire flux [kgC m-2 day-1]
 * m_litr4c_to_fire         double      fire flux [kgC m-2 day-1]
 * m_cwdc_to_fire           double      fire flux [kgC m-2 day-1]
 * leafc_transfer_to_leafc  double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * frootc_transfer_to_frootc
 *                          double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * livestemc_transfer_to_livestemc
 *                          double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * deadstemc_transfer_to_deadstemc
 *                          double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * livecrootc_transfer_to_livecrootc
 *                          double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * deadcrootc_transfer_to_deadcrootc
 *                          double      phenology flux from transfer pool
 *                                        [kgC m-2 day-1]
 * leafc_to_litr1c          double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * leafc_to_litr2c          double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * leafc_to_litr3c          double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * leafc_to_litr4c          double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * frootc_to_litr1c         double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * frootc_to_litr2c         double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * frootc_to_litr3c         double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * frootc_to_litr4c         double      leaf and fine root litterfall
 *                                        [kgC m-2 day-1]
 * leaf_day_mr              double      maintenance respiration flux
 *                                        [kgC m-2 day-1]
 * leaf_night_mr            double      maintenance respiration flux
 *                                        [kgC m-2 day-1]
 * froot_mr                 double      maintenance respiration flux
 *                                        [kgC m-2 day-1]
 * livestem_mr              double      maintenance respiration flux
 *                                        [kgC m-2 day-1]
 * livecroot_mr             double      maintenance respiration flux
 *                                        [kgC m-2 day-1]
 * psnsun_to_cpool          double      photosynthesis flux [kgC m-2 day-1]
 * psnshade_to_cpool        double      photosynthesis flux [kgC m-2 day-1]
 * cwdc_to_litr2c           double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * cwdc_to_litr3c           double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * cwdc_to_litr4c           double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr1_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr1c_to_soil1c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr2_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr2c_to_soil2c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr3c_to_litr2c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr4_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * litr4c_to_soil3c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil1_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil1c_to_soil2c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil2_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil2c_to_soil3c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil3_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil3c_to_soil4c         double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * soil4_hr                 double      litter decomposition flux
 *                                        [kgC m-2 day-1]
 * cpool_to_leafc           double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_leafc_storage   double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_frootc          double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_frootc_storage  double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_livestemc       double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_livestemc_storage
 *                          double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_deadstemc       double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_deadstemc_storage
 *                          double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_livecrootc      double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_livecrootc_storage
 *                          double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_deadcrootc      double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_deadcrootc_storage
 *                          double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_to_gresp_storage   double      daily allocation flux from current GPP
 *                                        [kgC m-2 day-1]
 * cpool_leaf_gr            double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_leaf_storage_gr    double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_leaf_gr         double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_froot_gr           double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_froot_storage_gr   double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_froot_gr        double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_livestem_gr        double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_livestem_storage_gr
 *                          double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_livestem_gr     double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_deadstem_gr        double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_deadstem_storage_gr
 *                          double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_deadstem_gr     double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_livecroot_gr       double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_livecroot_storage_gr
 *                          double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_livecroot_gr    double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_deadcroot_gr       double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * cpool_deadcroot_storage_gr
 *                          double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * transfer_deadcroot_gr    double      daily growth respiration flux
 *                                        [kgC m-2 day-1]
 * leafc_storage_to_leafc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * frootc_storage_to_frootc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * livestemc_storage_to_livestemc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * deadstemc_storage_to_deadstemc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * livecrootc_storage_to_livecrootc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * deadcrootc_storage_to_deadcrootc_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * gresp_storage_to_gresp_transfer
 *                          double      annual turnover of storage to transfer
 *                                        pools [kgC m-2 day-1]
 * livestemc_to_deadstemc   double      turnover of live wood to dead wood
 *                                        [kgC m-2 day-1]
 * livecrootc_to_deadcrootc double      turnover of live wood to dead wood
 *                                        [kgC m-2 day-1]
 ****************************************************************************/
typedef struct cflux_struct
{
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
    double          leafc_transfer_to_leafc;
    double          frootc_transfer_to_frootc;
    double          livestemc_transfer_to_livestemc;
    double          deadstemc_transfer_to_deadstemc;
    double          livecrootc_transfer_to_livecrootc;
    double          deadcrootc_transfer_to_deadcrootc;
    double          leafc_to_litr1c;
    double          leafc_to_litr2c;
    double          leafc_to_litr3c;
    double          leafc_to_litr4c;
    double          frootc_to_litr1c;
    double          frootc_to_litr2c;
    double          frootc_to_litr3c;
    double          frootc_to_litr4c;
    double          leaf_day_mr;
    double          leaf_night_mr;
    double          froot_mr;
    double          livestem_mr;
    double          livecroot_mr;
    double          psnsun_to_cpool;
    double          psnshade_to_cpool;
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
    double          leafc_storage_to_leafc_transfer;
    double          frootc_storage_to_frootc_transfer;
    double          livestemc_storage_to_livestemc_transfer;
    double          deadstemc_storage_to_deadstemc_transfer;
    double          livecrootc_storage_to_livecrootc_transfer;
    double          deadcrootc_storage_to_deadcrootc_transfer;
    double          gresp_storage_to_gresp_transfer;
    double          livestemc_to_deadstemc;
    double          livecrootc_to_deadcrootc;
} cflux_struct;


/*****************************************************************************
 * Nitrogen state variables (including sums for sources and sinks)
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * leafn                    double      leaf N [kgN m-2]
 * leafn_storage            double      leaf N [kgN m-2]
 * leafn_transfer           double      leaf N [kgN m-2]
 * frootn                   double      fine root N [kgN m-2]
 * frootn_storage           double      fine root N [kgN m-2]
 * frootn_transfer          double      fine root N [kgN m-2]
 * livestemn                double      live stem N [kgN m-2]
 * livestemn_storage        double      live stem N [kgN m-2]
 * livestemn_transfer       double      live stem N [kgN m-2]
 * deadstemn                double      dead stem N [kgN m-2]
 * deadstemn_storage        double      dead stem N [kgN m-2]
 * deadstemn_transfer       double      dead stem N [kgN m-2]
 * livecrootn               double      live coarse root N [kgN m-2]
 * livecrootn_storage       double      live coarse root N [kgN m-2]
 * livecrootn_transfer      double      live coarse root N [kgN m-2]
 * deadcrootn               double      dead coarse root N [kgN m-2]
 * deadcrootn_storage       double      dead coarse root N [kgN m-2]
 * deadcrootn_transfer      double      dead coarse root N [kgN m-2]
 * cwdn                     double      coarse woody debris N [kgN m-2]
 * litr1n                   double      litter labile N [kgN m-2]
 * litr2n                   double      litter unshielded cellulose N
 *                                        [kgN m-2]
 * litr3n                   double      litter shielded cellulose N [kgN m-2]
 * litr4n                   double      litter lignin N [kgN m-2]
 * soil1n                   double      microbial recycling pool N (fast)
 *                                        [kgN m-2]
 * soil2n                   double      microbial recycling pool N (medium)
 *                                        [kgN m-2]
 * soil3n                   double      microbial recycling pool N (slow)
 *                                        [kgN m-2]
 * soil4n                   double      recalcitrant SOM N (humus, slowest)
 *                                        [kgN m-2]
 * sminn                    double      soil mineral N [kgN m-2]
 * retransn                 double      plant pool of retranslocated N
 *                                        [kgN m-2]
 * npool                    double      temporary plant N pool [kgN m-2]
 * nfix_src                 double      SUM of biological N fixation [kgN m-2]
 * ndep_src                 double      SUM of N deposition inputs [kgN m-2]
 * nleached_snk             double      SUM of N leached [kgN m-2]
 * nvol_snk                 double      SUM of N lost to volatilization
 *                                        [kgN m-2]
 * fire_snk                 double      SUM of N lost to fire [kgN m-2]
 ****************************************************************************/
typedef struct nstate_struct
{
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
    double          nfix_src;
    double          ndep_src;
    double          nleached_snk;
    double          nvol_snk;
    double          fire_snk;
} nstate_struct;


/*****************************************************************************
 * Daily nitrogen flux variables
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * m_leafn_to_litr1n        double      mortality flux (kgN/m2/d)
 * m_leafn_to_litr2n        double      mortality flux (kgN/m2/d)
 * m_leafn_to_litr3n        double      mortality flux (kgN/m2/d)
 * m_leafn_to_litr4n        double      mortality flux (kgN/m2/d)
 * m_frootn_to_litr1n       double      mortality flux (kgN/m2/d)
 * m_frootn_to_litr2n       double      mortality flux (kgN/m2/d)
 * m_frootn_to_litr3n       double      mortality flux (kgN/m2/d)
 * m_frootn_to_litr4n       double      mortality flux (kgN/m2/d)
 * m_leafn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_frootn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_livestemn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_deadstemn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_livecrootn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_deadcrootn_storage_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_leafn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_frootn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_livestemn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_deadstemn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_livecrootn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_deadcrootn_transfer_to_litr1n
 *                          double      mortality flux (kgN/m2/d)
 * m_livestemn_to_litr1n    double      mortality flux (kgN/m2/d)
 * m_livestemn_to_cwdn      double      mortality flux (kgN/m2/d)
 * m_deadstemn_to_cwdn      double      mortality flux (kgN/m2/d)
 * m_livecrootn_to_litr1n   double      mortality flux (kgN/m2/d)
 * m_livecrootn_to_cwdn     double      mortality flux (kgN/m2/d)
 * m_deadcrootn_to_cwdn     double      mortality flux (kgN/m2/d)
 * m_retransn_to_litr1n     double      mortality flux (kgN/m2/d)
 * m_leafn_to_fire          double      fire flux (kgN/m2/d)
 * m_frootn_to_fire         double      fire flux (kgN/m2/d)
 * m_leafn_storage_to_fire  double      fire flux (kgN/m2/d)
 * m_frootn_storage_to_fire double      fire flux (kgN/m2/d)
 * m_livestemn_storage_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_deadstemn_storage_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_livecrootn_storage_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_deadcrootn_storage_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_leafn_transfer_to_fire double      fire flux (kgN/m2/d)
 * m_frootn_transfer_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_livestemn_transfer_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_deadstemn_transfer_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_livecrootn_transfer_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_deadcrootn_transfer_to_fire
 *                          double      fire flux (kgN/m2/d)
 * m_livestemn_to_fire      double      fire flux (kgN/m2/d)
 * m_deadstemn_to_fire      double      fire flux (kgN/m2/d)
 * m_livecrootn_to_fire     double      fire flux (kgN/m2/d)
 * m_deadcrootn_to_fire     double      fire flux (kgN/m2/d)
 * m_retransn_to_fire       double      fire flux (kgN/m2/d)
 * m_litr1n_to_fire         double      fire flux (kgN/m2/d)
 * m_litr2n_to_fire         double      fire flux (kgN/m2/d)
 * m_litr3n_to_fire         double      fire flux (kgN/m2/d)
 * m_litr4n_to_fire         double      fire flux (kgN/m2/d)
 * m_cwdn_to_fire           double      fire flux (kgN/m2/d)
 * leafn_transfer_to_leafn  double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * frootn_transfer_to_frootn
 *                          double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * livestemn_transfer_to_livestemn
 *                          double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * deadstemn_transfer_to_deadstemn
 *                          double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * livecrootn_transfer_to_livecrootn
 *                          double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * deadcrootn_transfer_to_deadcrootn
 *                          double      phenology flux from transfer pool
 *                                        (kgN/m2/d)
 * leafn_to_litr1n          double      litterfall flux (kgN/m2/d)
 * leafn_to_litr2n          double      litterfall flux (kgN/m2/d)
 * leafn_to_litr3n          double      litterfall flux (kgN/m2/d)
 * leafn_to_litr4n          double      litterfall flux (kgN/m2/d)
 * leafn_to_retransn        double      litterfall flux (kgN/m2/d)
 * frootn_to_litr1n         double      litterfall flux (kgN/m2/d)
 * frootn_to_litr2n         double      litterfall flux (kgN/m2/d)
 * frootn_to_litr3n         double      litterfall flux (kgN/m2/d)
 * frootn_to_litr4n         double      litterfall flux (kgN/m2/d)
 * ndep_to_sminn            double      deposition flux (kgN/m2/d)
 * nfix_to_sminn            double      deposition flux (kgN/m2/d)
 * cwdn_to_litr2n           double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * cwdn_to_litr3n           double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * cwdn_to_litr4n           double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * litr1n_to_soil1n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil1n_l1       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * litr2n_to_soil2n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil2n_l2       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * litr3n_to_litr2n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * litr4n_to_soil3n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil3n_l4       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * soil1n_to_soil2n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil2n_s1       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * soil2n_to_soil3n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil3n_s2       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * soil3n_to_soil4n         double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_soil4n_s3       double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * soil4n_to_sminn          double      litter and soil decomposition flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_l1s1       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_l2s2       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_l4s3       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_s1s2       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_s2s3       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_s3s4       double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_nvol_s4         double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_to_denitrif        double      denitrification (volatilization) flux
 *                                        (kgN/m2/d)
 * sminn_leached            double      leaching flux (kgN/m2/d)
 * retransn_to_npool        double      daily allocation flux (kgN/m2/d)
 * sminn_to_npool           double      daily allocation flux (kgN/m2/d)
 * npool_to_leafn           double      daily allocation flux (kgN/m2/d)
 * npool_to_leafn_storage   double      daily allocation flux (kgN/m2/d)
 * npool_to_frootn          double      daily allocation flux (kgN/m2/d)
 * npool_to_frootn_storage  double      daily allocation flux (kgN/m2/d)
 * npool_to_livestemn       double      daily allocation flux (kgN/m2/d)
 * npool_to_livestemn_storage
 *                          double      daily allocation flux (kgN/m2/d)
 * npool_to_deadstemn       double      daily allocation flux (kgN/m2/d)
 * npool_to_deadstemn_storage
 *                          double      daily allocation flux (kgN/m2/d)
 * npool_to_livecrootn      double      daily allocation flux (kgN/m2/d)
 * npool_to_livecrootn_storage
 *                          double      daily allocation flux (kgN/m2/d)
 * npool_to_deadcrootn      double      daily allocation flux (kgN/m2/d)
 * npool_to_deadcrootn_storage
 *                          double      daily allocation flux (kgN/m2/d)
 * leafn_storage_to_leafn_transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * frootn_storage_to_frootn_transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * livestemn_storage_to_livestemn_transfannual turnover of storage to transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * deadstemn_storage_to_deadstemn_transfannual turnover of storage to transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * livecrootn_storage_to_livecrootn_tranannual turnover of storage to transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * deadcrootn_storage_to_deadcrootn_tranannual turnover of storage to transfer
 *                          double      annual turnover of storage to transfer
 *                                        (kgN/m2/d)
 * livestemn_to_deadstemn   double      turnover of live wood to dead wood,
 *                                        with retranslocation (kgN/m2/d)
 * livestemn_to_retransn    double      turnover of live wood to dead wood,
 *                                        with retranslocation (kgN/m2/d)
 * livecrootn_to_deadcrootn double      turnover of live wood to dead wood,
 *                                        with retranslocation (kgN/m2/d)
 * livecrootn_to_retransn   double      turnover of live wood to dead wood,
 *                                        with retranslocation (kgN/m2/d)
 ****************************************************************************/
typedef struct nflux_struct
{
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
    double          leafn_transfer_to_leafn;
    double          frootn_transfer_to_frootn;
    double          livestemn_transfer_to_livestemn;
    double          deadstemn_transfer_to_deadstemn;
    double          livecrootn_transfer_to_livecrootn;
    double          deadcrootn_transfer_to_deadcrootn;
    double          leafn_to_litr1n;
    double          leafn_to_litr2n;
    double          leafn_to_litr3n;
    double          leafn_to_litr4n;
    double          leafn_to_retransn;
    double          frootn_to_litr1n;
    double          frootn_to_litr2n;
    double          frootn_to_litr3n;
    double          frootn_to_litr4n;
    double          ndep_to_sminn;
    double          nfix_to_sminn;
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
    double          sminn_to_nvol_l1s1;
    double          sminn_to_nvol_l2s2;
    double          sminn_to_nvol_l4s3;
    double          sminn_to_nvol_s1s2;
    double          sminn_to_nvol_s2s3;
    double          sminn_to_nvol_s3s4;
    double          sminn_to_nvol_s4;
    double          sminn_to_denitrif;
    double          sminn_leached;
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
    double          leafn_storage_to_leafn_transfer;
    double          frootn_storage_to_frootn_transfer;
    double          livestemn_storage_to_livestemn_transfer;
    double          deadstemn_storage_to_deadstemn_transfer;
    double          livecrootn_storage_to_livecrootn_transfer;
    double          deadcrootn_storage_to_deadcrootn_transfer;
    double          livestemn_to_deadstemn;
    double          livestemn_to_retransn;
    double          livecrootn_to_deadcrootn;
    double          livecrootn_to_retransn;
} nflux_struct;

/*****************************************************************************
 * Temporary nitrogen variables for reconciliation of decomposition
 * immobilization fluxes and plant growth N demands
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * mineralized              double      N mineralization [kgN m-2 day-1]
 * potential_immob          double      potential N immobilization [kgN m-2]
 * plitr1c_loss             double      potential loss from litter labile pool
 *                                        [kgN m-2 day-1]
 * pmnf_l1s1                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * plitr2c_loss             double      potential loss from litter unshielded
 *                                        pool [kgN m-2 day-1]
 * pmnf_l2s2                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * plitr4c_loss             double      potential loss from litter lignin pool
 *                                        [kgN m-2 day-1]
 * pmnf_l4s3                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * psoil1c_loss             double      potential loss from fast soil pool
 *                                        [kgN m-2 day-1]
 * pmnf_s1s2                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * psoil2c_loss             double      potential loss from medium soil pool
 *                                        [kgN m-2 day-1]
 * pmnf_s2s3                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * psoil3c_loss             double      potential loss from slow soil pool
 *                                        [kgN m-2 day-1]
 * pmnf_s3s4                double      potential mineral N flux
 *                                        [kgN m-2 day-1]
 * psoil4c_loss             double      potential loss from slowest soil pool
 *                                        [kgN m-2 day-1]
 * kl4                      double      decomposition rate of lignin litter
 *                                        pool [day -1]
 ****************************************************************************/
typedef struct ntemp_struct
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


/*****************************************************************************
 * Daily phenological data array
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * remdays_curgrowth        double      days left in current growth season
 * remdays_transfer         double      number of transfer days remaining
 * remdays_litfall          double      number of litfall days remaining
 * predays_transfer         double      number of transfer days previous
 * predays_litfall          double      number of litfall days previous
 ****************************************************************************/
typedef struct phenology_struct
{
    double          remdays_curgrowth;
    double          remdays_transfer;
    double          remdays_litfall;
    double          predays_transfer;
    double          predays_litfall;
} phenology_struct;


/*****************************************************************************
 * Ecophysiological variables
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * day_leafc_litfall_increment
 *                          double      rate leaf litfall [kgC m-2 day-1]
 * day_frootc_litfall_increment
 *                          double      rate froot litfall [kgC m-2 day-1]
 * day_livestemc_turnover_increment
 *                          double      rate livestem turnover [kgC m-2 day-1]
 * day_livecrootc_turnover_increment
 *                          double      rate livecroot turnover [kgC m-2 day-1]
 * annmax_leafc             double      annual maximum daily leaf C [kgC m-2]
 * annmax_frootc            double      annual maximum daily froot C [kgC m-2]
 * annmax_livestemc         double      annual maximum daily livestem C
 *                                        [kgC m-2]
 * annmax_livecrootc        double      annual maximum daily livecroot C
 *                                        [kgC m-2]
 * dsr                      double      number of days since rain
 * sun_proj_sla             double      sunlit projected SLA [m2 kgC-1]
 * shade_proj_sla           double      shaded projected SLA [m2 kgC-1]
 * psi                      double      water potential of soil and leaves [MPa]
 * dlmr_area_sun            double      sunlit leaf MR
 *                                        [umolC/m2 projected leaf area s-1]
 * dlmr_area_shade          double      shaded leaf MR
 *                                        [umolC/m2 projected leaf area s-1]
 * gl_t_wv_sun              double      leaf-scale conductance to transpired
 *                                        water [m s-1]
 * gl_t_wv_shade            double      leaf-scale conductance to transpired
 *                                        water [m s-1]
 * assim_sun                double      sunlit assimilation per unit pLAI
 *                                        [umol m-2 s-1]
 * assim_shade              double      shaded assimilation per unit pLAI
 *                                        [umol m-2 s-1]
 * t_scalar                 double      decomp temperature scalar [-]
 * w_scalar                 double      decomp water scalar [-]
 * rate_scalar              double      decomp combined scalar [-]
 * daily_gross_nmin         double      daily gross N mineralization
 *                                        [kgN m-2 d-1]
 * daily_gross_nimmob       double      daily gross N immobilization
 *                                        [kgN m-2 d-1]
 * daily_net_nmin           double      daily net N mineralization
 *                                        [kgN m-2 d-1]
 * fpi                      double      fraction of potential immobilization
 *                                        [-]
 * m_tmin                   double      freezing night temperature multiplier
 *                                        [-]
 * m_psi                    double      water potential multiplier [-]
 * m_co2                    double      atmospheric [CO2] multiplier [-]
 * m_ppfd_sun               double      PAR flux density multiplier [-]
 * m_ppfd_shade             double      PAR flux density multiplier [-]
 * m_vpd                    double      vapor pressure deficit multiplier [-]
 * m_final_sun              double      product of all other multipliers [-]
 * m_final_shade            double      product of all other multipliers [-]
 * ytd_maxplai              double      year-to-date maximum projected LAI [-]
 * dormant_flag             double      dormancy flag
 * days_active              double      number of days since last dormancy
 * onset_flag               double      onset flag
 * onset_counter            double      onset days counter
 * onset_gddflag            double      onset flag for growing degree day sum
 * onset_fdd                double      onset freezing degree days counter
 * onset_gdd                double      onset growing degree days
 * onset_swi                double      onset soil water index
 * offset_flag              double      offset flag
 * offset_counter           double      offset days counter
 * offset_fdd               double      offset freezing degree days counter
 * offset_swi               double      offset soil water index
 * lgsf                     double      long growing season factor (0-1)
 * bglfr                    double      background litterfall rate [s-1]
 * bgtr                     double      background transfer growth rate [s-1]
 * annavg_t2m               double      annual average 2m air temperature [K]
 * gpp                      double      GPP flux before downregulation
 *                                        [gC m-2 s-1]
 * prev_leafc_to_litter     double      previous timestep leaf C litterfall
 *                                        flux [gC m-2 s-1]
 * prev_frootc_to_litter    double      previous timestep froot C litterfall
 *                                        flux [gC m-2 s-1]
 * old_c_balance            double      previous timestep C balance
 *                                        [kgC m-2 day-1]
 * old_n_balance            double      previous timestep N balance
 *                                        [kgN m-2 day-1]
 ****************************************************************************/
typedef struct epvar_struct
{
    double          day_leafc_litfall_increment;
    double          day_frootc_litfall_increment;
    double          day_livestemc_turnover_increment;
    double          day_livecrootc_turnover_increment;
    double          annmax_leafc;
    double          annmax_frootc;
    double          annmax_livestemc;
    double          annmax_livecrootc;
    double          dsr;
    double          sun_proj_sla;
    double          shade_proj_sla;
    double          psi;
    double          dlmr_area_sun;
    double          dlmr_area_shade;
    double          gl_t_wv_sun;
    double          gl_t_wv_shade;
    double          assim_sun;
    double          assim_shade;
    double          t_scalar;
    double          w_scalar;
    double          rate_scalar;
    double          daily_gross_nmin;
    double          daily_gross_nimmob;
    double          daily_net_nmin;
    double          fpi;
    double          m_tmin;
    double          m_psi;
    double          m_co2;
    double          m_ppfd_sun;
    double          m_ppfd_shade;
    double          m_vpd;
    double          m_final_sun;
    double          m_final_shade;
    double          ytd_maxplai;
    double          dormant_flag;
    double          days_active;
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
    double          lgsf;
    double          bglfr;
    double          bgtr;
    double          annavg_t2m;
    double          gpp;
    double          prev_leafc_to_litter;
    double          prev_frootc_to_litter;
    double          old_c_balance;
    double          old_n_balance;
} epvar_struct;


/*****************************************************************************
 * Carbon state initialization structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * max_leafc                double      first-year displayed + stored leafc
 *                                        [kgC m-2]
 * max_stemc                double      first-year total stem carbon [kgC m-2]
 ****************************************************************************/
typedef struct cinit_struct
{
    double          max_leafc;
    double          max_stemc;
} cinit_struct;


/*****************************************************************************
 * Structure for the photosynthesis routine
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * c3                       int         set to 1 for C3 model, 0 for C4 model
 * pa                       double      atmospheric pressure [Pa]
 * co2                      double      atmospheric CO2 concentration [ppm]
 * t                        double      temperature [deg C]
 * lnc                      double      leaf N per unit sunlit leaf area
 *                                        [kg Nleaf m-2]
 * flnr                     double      fract. of leaf N in Rubisco
 *                                        [kg NRub/kg Nleaf]
 * ppfd                     double      PAR flux per unit sunlit leaf area 
 *                                        [umol m-2 s-1]
 * g                        double      conductance to CO2 [umol m-2 s-1 Pa-1]
 * dlmr                     double      day leaf maintenance respiration,
 *                                       projected area basis [umol m-2 s-1]
 * Ci                       double      intercellular CO2 concentration [Pa]
 * O2                       double      atmospheric O2 concentration [Pa]
 * Ca                       double      atmospheric CO2 concentration [Pa]
 * gamma                    double      CO2 compensation point, no Rd [Pa]
 * Kc                       double      MM constant carboxylation [Pa]
 * Ko                       double      MM constant oxygenation [Pa]
 * Vmax                     double      max rate carboxylation [umol m-2 s-1]
 * Jmax                     double      max rate electron transport
 *                                        [umol m-2 s-1]
 * J                        double      rate of RuBP regeneration
 *                                        [umol m-2 s-1]
 * Av                       double      carboxylation limited assimilation
 *                                        [umol m-2 s-1]
 * Aj                       double      RuBP regen limited assimilation
 *                                        [umol m-2 s-1]
 * A                        double      final assimilation rate [umol m-2 s-1]
 ****************************************************************************/
typedef struct psn_struct
{
    int             c3;
    double          pa;
    double          co2;
    double          t;
    double          lnc;
    double          flnr;
    double          ppfd;
    double          g;
    double          dlmr;
    double          Ci;
    double          O2;
    double          Ca;
    double          gamma;
    double          Kc;
    double          Ko;
    double          Vmax;
    double          Jmax;
    double          J;
    double          Av;
    double          Aj;
    double          A;
} psn_struct;

/*****************************************************************************
 * BGC summary structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * daily_npp                double      = GPP - Rmaint - Rgrowth
 *                                        [kgC m-2 day-1]
 * daily_nep                double      = NPP - Rheterotroph [kgC m-2 day-1]
 * daily_nee                double      = NEP - fire losses [kgC m-2 day-1]
 * daily_gpp                double      gross PSN source [kgC m-2 day-1]
 * daily_mr                 double      maintenance respiration
 *                                        [kgC m-2 day-1]
 * daily_gr                 double      growth respiration [kgC m-2 day-1]
 * daily_hr                 double      heterotrophic respiration
 *                                        [kgC m-2 day-1]
 * daily_fire               double      fire losses [kgC m-2 day-1]
 * daily_litfallc           double      total litterfall [kgC m-2 day-1]
 * cum_npp                  double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_nep                  double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_nee                  double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_gpp                  double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_mr                   double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_gr                   double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_hr                   double      Summed over entire simulation
 *                                        [kgC m-2]
 * cum_fire                 double      Summed over entire simulation
 *                                        [kgC m-2]
 * vegc                     double      total vegetation C [kgC m-2]
 * litrc                    double      total litter C [kgC m-2]
 * soilc                    double      total soil C [kgC m-2]
 * totalc                   double      total of vegc, litrc, and soilc
 *                                        [kgC m-2]
 ****************************************************************************/
typedef struct summary_struct
{
    double          daily_npp;
    double          daily_nep;
    double          daily_nee;
    double          daily_gpp;
    double          daily_mr;
    double          daily_gr;
    double          daily_hr;
    double          daily_fire;
    double          daily_litfallc;
    double          cum_npp;
    double          cum_nep;
    double          cum_nee;
    double          cum_gpp;
    double          cum_mr;
    double          cum_gr;
    double          cum_hr;
    double          cum_fire;
    double          vegc;
    double          litrc;
    double          soilc;
    double          totalc;
} summary_struct;
#endif

/*****************************************************************************
 * Element structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * node                     int[]       nodes of triagular element
 *                                        (Counterclock-wise)
 * nabr                     int[]       neighbor elements (neighbor i shares
 *                                        edge i (0: on boundary)
 * ind                      int         element index
 ****************************************************************************/
typedef struct elem_struct
{
    int             node[NUM_EDGE];     /* Counterclock-wise */
    int             nabr[NUM_EDGE];     /* neighbor i shares edge i
                                         * (0: on boundary) */
    int             ind;

    attrib_struct   attrib;
    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
    epconst_struct  epc;
    ic_struct       ic;
    bc_struct       bc;
    wstate_struct   ws;
    wstate_struct   ws0;
    wflux_struct    wf;
    estate_struct   es;
    eflux_struct    ef;
#ifdef _NOAH_
    wflux_struct    avgwf;
#endif
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
#endif
#ifdef _BGC_
    bgcic_struct    restart_input;
    bgcic_struct    restart_output;
    stor_struct     stor;
    cinit_struct    cinit;
    cstate_struct   cs;
    cflux_struct    cf;
    nstate_struct   ns;
    nflux_struct    nf;
    psn_struct      psn_sun;
    psn_struct      psn_shade;
    ntemp_struct    nt;
    summary_struct  summary;
    phenology_struct phen;
    epvar_struct    epv;
#endif
} elem_struct;
#endif
