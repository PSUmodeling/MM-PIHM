#include "Cycles.h"

void NitrogenTransformation (int y, int doy, SoilStruct *Soil, const CommunityStruct *Community, const ResidueStruct *Residue, const WeatherStruct *Weather, const SoilCarbonStruct *SoilCarbon)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i                    int
     * Profile_N_Nitrified  double      N daily nitrification [Mg N/ha]
     * Profile_N2O_Nitri    double      N daily N2O from nitrification [Mg/ha]
     * Profile_N_Denitrified
     *                      double      N daily denitrification [Mg N/ha]
     * Profile_N2O_Denit    double      N daily denitrification in the form of
     *                                    N2O [Mg N/ha]
     * Profile_NH4_Volatilization
     *                      double      NH4 volatilization [Mg N/ha]
     */
    int             i;
    double          Profile_N_Nitrified;
    double          Profile_N2O_Nitri;
    double          Profile_N_Denitrified;
    double          Profile_N2O_Denit;
    double          Profile_NH4_Volatilization;

    Profile_N_Nitrified = 0.0;
    Profile_N2O_Nitri = 0.0;
    Profile_N_Denitrified = 0.0;
    Profile_N2O_Denit = 0.0;
    Profile_NH4_Volatilization = 0.0;

    Soil->NH4_Nitrification = 0.0;
    Soil->NO3_Denitrification = 0.0;
    Soil->N2O_Denitrification = 0.0;
    Soil->NH4_Volatilization = 0.0;

    Denitrification (&Profile_N_Denitrified, &Profile_N2O_Denit, Soil, SoilCarbon);
    Nitrification (&Profile_N_Nitrified, &Profile_N2O_Nitri, Soil, SoilCarbon);
    Volatilization (y, doy, &Profile_NH4_Volatilization, Soil, Community, Residue, Weather);

    Soil->NO3Profile = 0.0;
    Soil->NH4Profile = 0.0;
    for (i = 0; i < Soil->totalLayers; i++)
    {
        Soil->NO3Profile += Soil->NO3[i];
        Soil->NH4Profile += Soil->NH4[i];
    }

    Soil->NH4_Nitrification = Profile_N_Nitrified;
    Soil->N2O_Nitrification = Profile_N2O_Nitri;
    Soil->NO3_Denitrification = Profile_N_Denitrified;
    Soil->N2O_Denitrification = Profile_N2O_Denit;
    Soil->NH4_Volatilization = Profile_NH4_Volatilization;
}

void Nitrification (double *Profile_N_Nitrified, double *Profile_N2O_Nitrified, SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i                    int
     * NH4_Nitrified        double      N daily nitrification [Mg/ha]
     * N2O_Nitrified        double      N daily N2O from denitrification and
     *                                    nitrification [Mg/ha]
     * pH_Factor            double      pH function controlling nitrification
     * NO3NH4Ratio          double      Ratio of NO3 to NH4
     * ratioFactor          double      Nitrification control by NO3/NH4 ratio
     * AirContent           double      Porosity - watercontent [m3/m3]
     * AirFactor            double      Nitrification control by air content
     * N2O_Fraction         double      Fraction of nitrification released as
     *                                    N2O
     * TempFactor           double      Nitrification control by temperature
     */
    int             i;
    double          NH4_Nitrified;
    double          N2O_Nitrified;
    double          pH_Factor = 1.0;
    double          NO3NH4Ratio;
    double          ratioFactor;
    double          AirContent;
    double          AirFactor;
    double          N2O_Fraction;
    double          TempFactor;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        if (Soil->NH4[i] > 0.0)
        {
            NO3NH4Ratio = Soil->NO3[i] / Soil->NH4[i];
            ratioFactor = 1.0 / (1.0 + pow (NO3NH4Ratio / 8.0, 6));
            AirContent = Soil->Porosity[i] - Soil->waterContent[i];
            AirFactor = 1.0 - 1.0 / (1.0 + pow (AirContent / 0.1, 3));
            N2O_Fraction = N2OFractionNitrification (AirContent);
            TempFactor = TemperatureFunction (Soil->soilTemperature[i]);
            NH4_Nitrified = Soil->NH4[i] * NITRIFICATION_CONSTANT * ratioFactor * pH_Factor * AirFactor * TempFactor;
            N2O_Nitrified = N2O_Fraction * NH4_Nitrified;
            Soil->NH4[i] -= (NH4_Nitrified + N2O_Nitrified);
            Soil->NO3[i] += NH4_Nitrified;
            *Profile_N_Nitrified = *Profile_N_Nitrified + NH4_Nitrified;
            *Profile_N2O_Nitrified = *Profile_N2O_Nitrified + N2O_Nitrified;

            Soil->n2o[i] = N2O_Nitrified;
        }
        else
        {
            NH4_Nitrified = 0.0;
            N2O_Nitrified = 0.0;
        }
    }
}

void Denitrification (double *Profile_N_Denitrified, double *Profile_N2O_Denitrified, SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * N_Denit              double      Nitrogen daily denitrification (N2 +
     *                                    N2O) [Mg N/ha]
     * N2O_Emission         double      Nitrogen daily denitrification as N2O
     *                                    [kg/m2]
     * N2O_Fraction         double      Ratio of N2O to total denitrification
     * Soil_Mass            double      [Mg/ha]
     * NO3_Conc             double      [Mg NO3 / Mg dry soil]
     * NO3_Factor           double      Nitrate concentration control of
     *                                    denitrification
     * Res_Factor           double      This respiration factor considers
     *                                    temperature control already
     * Oxy_Factor           double      Oxygen control of denitrification,
     *                                    using porosity occupied by air as a
     *                                    surrogate
     * rr1                  double      Carbon respired [Mg C / Mg dry soil]
     * AirVol               double      Fractional volume of air in the soil
     *                                    [m3/m3]
     * cc1                  double      Compute coefficient for
     *                                    denitrification based on clay
     *                                    concentration
     * cc2                  double      Coefficient of denitrification curve
     *                                    response to aereation
     * i                    int         Loop counter
     */

    double          N_Denit;
    double          N2O_Emission;
    double          N2O_Fraction;
    double          Soil_Mass;
    double          NO3_Conc;
    double          NO3_Factor;
    double          Res_Factor;
    double          Oxy_Factor;
    double          rr1;
    double          AirVol;
    double          cc1;
    const double    cc2 = 60.0;
    int             i;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        N_Denit = 0.0;
        N2O_Emission = 0.0;
        AirVol = Soil->Porosity[i] - Soil->waterContent[i];

        if (Soil->NO3[i] > 1e-6 && AirVol < 0.25)
        {
            cc1 = 0.9 + 0.1 * Soil->Clay[i];
            Oxy_Factor = 1.0 / (1.0 + pow ((1.0 - AirVol) / cc1, -cc2));

            Soil_Mass = Soil->BD[i] * Soil->layerThickness[i] * 10000.0;    /* converted to Mg soil/ha */
            rr1 = SoilCarbon->carbonRespired[i] / Soil_Mass;
            Res_Factor = rr1 / 0.00005 < 1.0 ? rr1 / 0.00005 : 1.0;

            NO3_Conc = Soil->NO3[i] / Soil_Mass;
            NO3_Factor = NO3_Conc / (NO3_Conc + DENITRIFICATION_HALF_RATE);

            //N_Denit = POTENTIAL_DENITRIFICATION * Soil_Mass * Oxy_Factor * Res_Factor * NO3_Factor
            N_Denit = POTENTIAL_DENITRIFICATION * Soil_Mass * pow (Oxy_Factor, 0.5) * Res_Factor * NO3_Factor;
            N_Denit = N_Denit < Soil->NO3[i] ? N_Denit : Soil->NO3[i];
            N2O_Fraction = NO3_Factor * (1.0 - pow (Oxy_Factor, 0.5)) * (1.0 - pow (Res_Factor, 0.33));
            N2O_Emission = N_Denit * N2O_Fraction;

            Soil->NO3[i] -= N_Denit;
        }

        *Profile_N_Denitrified += N_Denit;
        *Profile_N2O_Denitrified += N2O_Emission;

        Soil->n2o[i] = N2O_Emission;
    }
}

void Volatilization (int y, int doy, double *Profile_NH4_Volatilization, SoilStruct *Soil, const CommunityStruct *Community, const ResidueStruct *Residue, const WeatherStruct *Weather)
{
    /*
     * This subroutine uses an empirical approach to estimate the amount of
     * NH4 that is subject to volatilization and applies Henry equilibrium to
     * estimate NH3 volatilization using turbulent conductance of the
     * atmosphere and a proxy for conductance from the soil surface to the
     * canopy exchange surface an Eulerian approach with daily time step can
     * cause excess volatilization in a day
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i                    int         Loop counter
     * layerTop             double[]    [m]
     * layerBottom          double[]    [m]
     * layerMidpoint        double      [m]
     * DepthFactor          double      Depth factor affection fraction of
     *                                    ammonium volatilizable (no diffusion
     *                                    processes in soil simulated)
     * CEC                  double      CEC, estimated from clay and soil
     *                                    organic carbon
     * CECFactor            double      CEC factor controlling fraction of
     *                                    ammonia that is available for
     *                                    volatilization
     * fVol                 double      Fraction of layer ammonium that can
     *                                    volatilize
     * NH4Volatilizable     double      NH4 of a layer subject to
     *                                    volatilization
     * NH4Volatilized       double      NH4 volatilized from a given layer
     * Tavg                 double      Weighted daily average temperature
     *                                    (favoring daytime temperature)
     * pAtm                 double      Atmospheric pressure [Pa]
     * AMD                  double      Air molar density [mol/m3]
     * GBL                  double      Molar atmospheric boundary layer
     *                                    conductance [mol/m3/s]
     * GG1                  double      NH3 flux correction based on flat
     *                                    residue cover, 0 - 1 temporary
     * GG2                  double      NH3 flux correction based on canopy
     *                                    cover,  0 - 1 temporary
     * GG3                  double      Product GBL*GG2*GG3
     * pK                   double
     * pH                   double
     * henrysCoeff          double
     * henrysConst          double      Henry's constant, [Litre Pa / mol]
     * waterVolume          double
     * NH4Conc              double
     * NH3Conc              double
     * NH3MolarFraction     double      Molar fraction
     * profile_volatilization
     *                      double
     */

    int             i = 0;
    double          layerTop[Soil->totalLayers];
    double          layerBottom[Soil->totalLayers];
    double          layerMidpoint;
    double          DepthFactor;
    double          CEC;
    double          CECFactor;
    double          fVol;
    double          NH4Volatilizable;
    double          NH4Volatilized;
    double          Tavg;
    double          pAtm;
    double          AMD;
    double          GBL;
    double          GG1;
    double          GG2;
    double          GG3;
    double          pK;
    double          pH;
    double          henrysCoeff;
    double          henrysConst;
    double          waterVolume;
    double          NH4Conc;
    double          NH3Conc;
    double          NH3MolarFraction;
    double          profile_volatilization;

    pH = 6.5;

    profile_volatilization = *Profile_NH4_Volatilization;

    Tavg = 273.15 + 0.67 * Weather->tMax[y][doy - 1] + 0.33 * Weather->tMin[y][doy - 1];
    pAtm = Weather->atmosphericPressure * 1000.0;
    AMD = AirMolarDensity (Tavg, pAtm);
    GBL = BoundaryLayerConductance (Community->svRadiationInterception, Residue->stanResidueMass, Weather->wind[y][doy - 1], AMD);
    GG1 = 1.0 - 0.85 * pow (Community->svRadiationInterception, 3.0);
    GG2 = 0.95 * pow (Residue->flatResidueTau, 2.0);
    GG3 = GBL * GG1 * GG2;

    pK = 0.09018 + 2729.92 / Tavg;
    henrysCoeff = pow (10.0, 1477.7 / Tavg - 1.69);
    henrysConst = 1000.0 * 8.3143 * Tavg / henrysCoeff; /* 1000 converts from
                                                         * m3 to liters */

    for (i = 0; i < Soil->totalLayers; i++)
    {
        if (i > 0)
        {
            layerBottom[i] = layerBottom[i - 1] + Soil->layerThickness[i];
            layerTop[i] = layerBottom[i - 1];
        }
        else
        {
            layerBottom[i] = Soil->layerThickness[i];
            layerTop[i] = 0.0;
        }

        CEC = 0.7 * Soil->Clay[i] * 100.0;
        CEC = CEC < 13.0 * Soil->SOC_Conc[i] * 100.0 ? CEC : 13.0 * Soil->SOC_Conc[i] * 100.0;
        CECFactor = 0.2 + 0.8 * exp (-0.08 * CEC);

        layerMidpoint = 0.5 * (layerTop[i] + layerBottom[i]);
        DepthFactor = 1.0 / 3.0 * (VolatilizationDepthFunction (layerTop[i]) + VolatilizationDepthFunction (layerMidpoint) + VolatilizationDepthFunction (layerBottom[i]));

        fVol = CECFactor * DepthFactor;
        NH4Volatilizable = fVol * Soil->NH4[i];

        waterVolume = Soil->waterContent[i] * Soil->layerThickness[i] * WATER_DENSITY * 10000.0;    /* m3/ha */
        NH4Conc = NH4Volatilizable * (18.0 / 14.0) / waterVolume;   /* Mg/m3; 18/14 converts mass of N to mass of NH4 */
        NH3Conc = NH4Conc / (1.0 + pow (10.0, pK - pH));    /* Mg/m3 */
        NH3MolarFraction = henrysConst * (NH3Conc / 0.000017) / pAtm;   /* 0.000017 = Mg/mol of NH3 */

        NH4Volatilized = GG3 * NH3MolarFraction * 86400.0 * 0.000017 * 10000.0 * (14.0 / 17.0); /* Mg NH3 / ha / day; 14/17 converts mass of N to mass of NH4 */

        Soil->NH4[i] -= NH4Volatilized;
        profile_volatilization += NH4Volatilized;
    }

    *Profile_NH4_Volatilization = profile_volatilization;
}

double N2OFractionNitrification (double air)
{
    /*
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * N2OFraction_Function double
     * Q                    double
     * f_Base               double
     * f_Max                double
     * a_Min                double
     * a_Opt                double
     * a_Max                double
     */

    double          N2OFraction_Function;
    double          Q;
    const double    f_Base = 0.0025;
    const double    f_Max = 0.1;
    const double    a_Min = 0.0;
    const double    a_Opt = 0.05;
    const double    a_Max = 0.1;

    if (air < 0.0 || air > a_Max)
        N2OFraction_Function = f_Base;
    else
    {
        Q = (a_Min - a_Opt) / (a_Opt - a_Max);
        N2OFraction_Function = f_Base + f_Max * (pow (air - a_Min, Q) * (a_Max - air)) / (pow (a_Opt - a_Min, Q) * (a_Max - a_Opt));
        if (N2OFraction_Function > 1.0)
            N2OFraction_Function = 1.0;
    }

    return (N2OFraction_Function);
}

double pHFunction (double pH)
{
    /*
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * pH_Min               double
     * pH_Max               double
     */

    double          pH_Min = 3.5;
    double          pH_Max = 6.5;

    return (pH - pH_Min) / (pH_Max - pH_Min);
}

double VolatilizationDepthFunction (double depth)
{
    return (1.0 / (1.0 + pow (depth / 0.03, 3.5)));
}

double AirMolarDensity (double T, double P)
{
    const double    Rgas = 8.3143;
    return P / (Rgas * T);
}

double BoundaryLayerConductance (double RI, double RM, double WS, double AMD)
{
    /* RI = fractional solar radiation interception, 0 - 1
     * RM = mass of standing residue, Mg/ha
     * AMD = air molar density, mol/m3
     * WS = wind speed at 2 meters, m/s */
    /*
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * CropHeight           double
     * SoilHeight           double
     * ResidueHeight        double
     * h                    double
     * d                    double      Zero plane displacement (at
     *                                    d' = d + Zm, wind speed profile
     *                                    extrapolates to zero)
     * Zm                   double      Momentum roughness length
     * Zs                   double      Scalar roughness length
     */

    double          CropHeight;
    double          SoilHeight;
    double          ResidueHeight;

    double          h;
    double          d;          /* zero plane displacement (at d' = d + Zm, wind speed profile extrapolates to zero) */
    double          Zm;         /* momentum roughness length */
    double          Zs;         /* scalar roughness length */

    CropHeight = RI / (RI + 0.5);   /* temporary, until a function for crop height is incorporated; assumed max = 1 m */
    SoilHeight = 0.05;          /* temporary, until something better is available */
    ResidueHeight = 0.3 * RM / (RM + 2.0);  /* 2 in Mg/ha or residue, assumed max = 0.3 m */

    h = CropHeight > SoilHeight ? CropHeight : SoilHeight;
    h = h > ResidueHeight ? h : ResidueHeight;
    d = 0.65 * h;
    Zm = 0.1 * h;
    Zs = 0.2 * Zm;

    return pow (0.41, 2.0) * AMD * WS / (log ((2.0 - d) / Zm) * log ((2.0 - d) / Zs));
}
