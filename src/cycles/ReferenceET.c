#include "Cycles.h"

double CalculatePMET (double lat, double pAtm, double screeningHeight, double Tmax, double Tmin, double sRad, double rhMax, double rhMin, double wind, double doy)
{
    /*
     * Calculate potential ET using Penman Monteith
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * PMET		    double
     * rc		    double	Surface resistance to vapor flow
     *					  (day/m) (0.00081 for 350-400 ppm)
     * CP		    double	Specific heat of dry air (J/kg/C)
     * r_gas		    double	Specific gas constant for dry air
     *					  (kJ/kg/degree C) (recall that Patm
     *					  is in kPa)
     * esTmax		    double	Saturation vapor pressure at Tmax (kPa)
     * esTmin		    double	Saturation vapor pressure at Tmin (kPa)
     * ea		    double	Actual air vapor pressure (kPa)
     * vpd		    double
     * potRad		    double
     * netRad		    double
     * Tave		    double
     * esTave		    double
     * delta		    double
     * gamma		    double
     * lambda		    double
     * Tkv		    double
     * volCP		    double
     * ra		    double
     * aeroTerm		    double
     * radTerm		    double
     */
    double          PMET;
    const double    rc = 0.00081;
    const double    CP = 0.001013;
    const double    r_gas = 0.28704;
    double          esTmax;
    double          esTmin;
    double          ea;
    double          vpd;
    double	    potRad;
    double	    netRad;
    double	    Tave;
    double	    esTave;
    double	    delta;
    double	    gamma;
    double	    lambda;
    double	    Tkv;
    double	    volCP;
    double	    ra;
    double	    aeroTerm;
    double	    radTerm;

    Tave = (Tmax + Tmin) / 2.0;
    esTave = SaturatedVaporPressure (Tave);
    esTmax = SaturatedVaporPressure (Tmax);
    esTmin = SaturatedVaporPressure (Tmin);
    ea = 0.5 * (esTmin * rhMax + esTmax * rhMin) / 100.0;
    vpd = (esTmax + esTmin) / 2.0 - ea;
    potRad = PotentialRadiation (lat, doy);
    netRad = NetRadiation (potRad, sRad, ea, Tmax, Tmin);

    /* Aerodynamic resistance to vapor flow (day/m) */
    ra = Aero_Res (wind, screeningHeight);

    /* Slope of saturated vapor pressure vs temperature function (kPa/C) */
    delta = 4098.0 * esTave / (pow (Tave + 237.3, 2.0));

    /* Latent heat of vaporization (MJ/kg) */
    lambda = 2.501 - 0.002361 * Tave;

    /* Psychrometric constant (kPaC) */
    gamma = CP * pAtm / (0.622 * lambda);

    /* Approximates virtual temperature (K) */
    Tkv = 1.01 * (Tave + 273.15);

    /* CP * AirDensity (J/kg * kg/m3) */
    volCP = CP * pAtm / (r_gas * Tkv);

    /* Aerodynamic term (MJ/m2) */
    aeroTerm = (volCP * vpd / ra) / (delta + gamma * (1.0 + rc / ra));

    /* Radiation term (MJ/m2) */
    radTerm = delta * netRad / (delta + gamma * (1.0 + rc / ra));

    /* Potential ET (kg water/m2 or mm) (water density = 1 Mg/m3) */
    PMET = (aeroTerm + radTerm) / lambda;

    /* Preventing a negative value usually small and indicative of
     * condensation */
    if (PMET < 0)
        PMET = 0.001;

    return (PMET);
}

double SaturatedVaporPressure (double T)
{
    return (0.6108 * exp (17.27 * T / (T + 237.3)));
}

double PotentialRadiation (double Lat, int doy)
{
    /*
     * Calculate potential radiation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * Solar_Constant	    double
     * Lat_Rad		    double
     * DR		    double
     * SolDec		    double
     * SunsetHourAngle	    double
     * Term		    double
     */
    const double    Solar_Constant = 118.08;
    double          Lat_Rad;
    double	    DR;
    double	    SolDec;
    double	    SunsetHourAngle;
    double	    Term;

    Lat_Rad = Lat * PI / 180.0;
    DR = 1.0 + 0.033 * cos (2.0 * PI * doy / 365.0);
    SolDec = 0.409 * sin (2.0 * PI * doy / 365.0 - 1.39);
    SunsetHourAngle = acos (-tan (Lat_Rad) * tan (SolDec));
    Term = SunsetHourAngle * sin (Lat_Rad) * sin (SolDec) + cos (Lat_Rad) * cos (SolDec) * sin (SunsetHourAngle);

    return (Solar_Constant * DR * Term / PI);
}

double NetRadiation (double Pot_Rad, double Solar_Rad, double Actual_VP, double TMax, double TMin)
{
    /*
     * Calculate net radiation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * Rns		    double	Shartwave net radiation
     * F_Cloud		    double	Cloud factor
     * F_Hum		    double	Humidity factor
     * LWR		    double	Isothermal longwave net radiation
     * Rnl		    double
     * Albedo		    double	Albedo (-)
     */
    double          Rns;
    double	    F_Cloud;
    double	    F_Hum;
    double	    LWR;
    double	    Rnl;
    const double    Albedo = 0.23;

    /* Calculate shortwave net radiation */
    Rns = (1.0 - Albedo) * Solar_Rad;

    /* Calculate cloud factor */
    F_Cloud = 1.35 * (Solar_Rad / (Pot_Rad * 0.75)) - 0.35;

    /* Calculate humidity factor */
    F_Hum = (0.34 - 0.14 * sqrt (Actual_VP));

    /* Calculate isothermal LW net radiation */
    LWR = 86400.0 * 5.67e-8 * 0.5 * (pow (TMax + 273.15, 4.0) + pow (TMin + 273.15, 4.0)) / 1.0e6; 

    Rnl = LWR * F_Cloud * F_Hum;

    return (Rns - Rnl);
}

double Aero_Res (double uz, double z)
{
    /*
     * Calculate aerodynamic resistance
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * u2		    double
     * d		    double	Zero plane displacement height (m)
     * zom		    double	Roughness length for momentum transfer
     *					  (m)
     * zoh		    double	Roughness length for heat and vapor
     *					  transfer (m)
     * zm		    double	Height of wind measurement (m)
     * zh		    double	Height of vapor measurement (m)
     * VK		    double	von Karman's constant
     * r_a		    double	Aerodynamic resistance
     */
    double          u2;
    double	    d;
    double	    zom;
    double	    zoh;
    double	    zm;
    double	    zh;
    const double    VK = 0.41;
    double          r_a;

    if (uz == 0.0)
        uz = 0.00001;

    if (z == 2.0)
        u2 = uz;
    else
        u2 = uz * (4.87 / (log (67.8 * z - 5.42)));

    u2 = u2 * 86400.0;          /* convert to m/day */
    d = 0.08;
    zom = 0.01476;
    zoh = 0.001476;
    zm = 2.0;
    zh = 2.0;

    r_a = log ((zm - d) / zom) * log ((zh - d) / zoh) / (VK * VK * u2);

    return (r_a);
}
