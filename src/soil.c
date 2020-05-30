#include "pihm.h"

double KrFunc(double beta, double satn)
{
    return sqrt(satn) *
        (1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta)) *
        (1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta));
}

double FieldCapacity(double beta, double kv, double smcmax,
    double smcmin)
{
    /*
     * Solve field capacity using Newton's method
     * Field capacity is defined as (Chen and Dudhia 2001 MWR)
     * Theta_ref = Theta_s * (1 / 3 + 2 / 3 * satn_ref)
     * where satn_ref is the soil saturation ratio when the soil hydraulic
     * conductivity reaches 5.70E-9 m/s
     */
    double          satn;
    double          satnk;
    double          denom;
    double          df;
    double          dsatn = 999.0;
    double          mx;
    double          ftheta;
    int             n;
    double          smcref;
    const double    KFC = 5.79E-9;
    const double    TOL = 1.0E-3;

    mx = 1.0 - 1.0 / beta;

    n = 0;

    satn = 0.75;

    while (n < 10 && dsatn > TOL)
    {
        n++;

        df = KrFunc(beta, satn) - KFC / kv;

        ftheta = 1.0 - pow(satn, 1.0 / mx);

        denom = 0.5 * pow(satn, -0.5) * pow(1.0 - pow(ftheta, mx), 2.0) +
            2.0 * pow(satn, 1.0 / mx - 0.5) *
            (pow(ftheta, mx - 1.0) - pow(ftheta, 2.0 * mx - 1.0));

        satnk = satn - df / denom;
        satnk = (satnk > 1.0 - 1.0E-3) ? 1.0 - 1.0E-3 : satnk;
        satnk = (satnk < SATMIN) ? SATMIN : satnk;

        dsatn = fabs(satnk - satn);

        satn = satnk;
    }

    if (dsatn > TOL)
    {
        satn = 0.75;
    }

    smcref = (1.0 / 3.0 + 2.0 / 3.0 * satn) * (smcmax - smcmin) + smcmin;

    return smcref;
}

double WiltingPoint(double smcmax, double smcmin, double alpha, double beta)
{
    const double    PSIW = 200.0;

    return 0.5 * (smcmax - smcmin) *
        pow(1.0 / (1.0 + pow(PSIW * alpha, beta)), 1.0 - 1.0 / beta) + smcmin;
}

int SoilTex(double silt, double clay)
{
    /*
     * Define soil texture using USDA Textural Classes
     */
    int             texture = -999;
    double          sand;

    silt /= 100.0;
    clay /= 100.0;

    sand = 1.0 - silt - clay;

    if (silt < 0.0 || silt > 1.0)
    {
        pihm_printf(VL_ERROR,
            "Error: Silt percentage (%lf) out of range.\n", silt * 100.0);
        pihm_printf(VL_ERROR, "Please check your soil input file.\n");
        pihm_exit(EXIT_FAILURE);
    }
    if (clay < 0.0 || clay > 1.0)
    {
        pihm_printf(VL_ERROR,
            "Error: Clay percentage (%lf) out of range.\n", clay * 100.0);
        pihm_printf(VL_ERROR, "Please check your soil input file.\n");
        pihm_exit(EXIT_FAILURE);
    }
    if (sand < 0.0 || sand > 1.0)
    {
        pihm_printf(VL_ERROR,
            "Error: Sand percentage (%lf) out of range.\n", sand * 100.0);
        pihm_printf(VL_ERROR, "Please check your soil input file.\n");
        pihm_exit(EXIT_FAILURE);
    }

    if (silt + 1.5 * clay < 0.15)
    {
        texture = SAND;
    }
    else if (silt + 1.5 * clay >= 0.15 && silt + 2.0 * clay < 0.3)
    {
        texture = LOAMY_SAND;
    }
    else if ((clay >= 0.07 && clay <= 0.2 && sand > 0.52 &&
            silt + 2.0 * clay >= 0.3) || (clay < 0.07 && silt < 0.5 &&
            silt + 2.0 * clay >= 0.3))
    {
        texture = SANDY_LOAM;
    }
    else if (clay >= 0.07 && clay <= 0.27 && silt >= 0.28 && silt < 0.5 &&
        sand <= 0.52)
    {
        texture = LOAM;
    }
    else if ((silt > 0.5 && clay > 0.12 && clay < 0.27) ||
        (silt >= 0.5 && silt < 0.8 && clay < 0.12))
    {
        texture = SILT_LOAM;
    }
    else if (silt >= 0.8 && clay < 0.12)
    {
        texture = SILT;
    }
    else if (clay >= 0.2 && clay < 0.35 && silt < 0.28 && sand > 0.45)
    {
        texture = SANDY_CLAY_LOAM;
    }
    else if (clay >= 0.27 && clay < 0.4 && sand > 0.2 && sand <= 0.45)
    {
        texture = CLAY_LOAM;
    }
    else if (clay >= 0.27 && clay < 0.4 && sand <= 0.2)
    {
        texture = SILTY_CLAY_LOAM;
    }
    else if (clay >= 0.35 && sand > 0.45)
    {
        texture = SANDY_CLAY;
    }
    else if (clay >= 0.4 && silt >= 0.4)
    {
        texture = SILTY_CLAY;
    }
    else if (clay >= 0.4 && sand <= 0.45 && silt < 0.4)
    {
        texture = CLAY;
    }
    else
    {
        pihm_printf(VL_ERROR, "Error: Soil texture %d not defined.\n", texture);
        pihm_exit(EXIT_FAILURE);
    }

    return texture;
}

double Qtz(int texture)
{
    double          qtz = 0.0;

    switch (texture)
    {
        case SAND:
            qtz = 0.92;
            break;
        case LOAMY_SAND:
            qtz = 0.82;
            break;
        case SANDY_LOAM:
            qtz = 0.60;
            break;
        case SILT_LOAM:
            qtz = 0.25;
            break;
        case SILT:
            qtz = 0.10;
            break;
        case LOAM:
            qtz = 0.40;
            break;
        case SANDY_CLAY_LOAM:
            qtz = 0.60;
            break;
        case SILTY_CLAY_LOAM:
            qtz = 0.10;
            break;
        case CLAY_LOAM:
            qtz = 0.35;
            break;
        case SANDY_CLAY:
            qtz = 0.52;
            break;
        case SILTY_CLAY:
            qtz = 0.10;
            break;
        case CLAY:
            qtz = 0.25;
            break;
        default:
            pihm_printf(VL_ERROR, "Error: Soil texture %d not defined.\n",
                texture);
            pihm_exit(EXIT_FAILURE);
            break;
    }

    return qtz;
}

/*
 * Pedotransfer functions to calculate soil hydraulic properties
 * Wosten et al. 1999 Geoderma Table 5
 */
double PtfKv(double silt, double clay, double om, double bd, int topsoil)
{
    double          kv;

    /* Calculate Kv in cm/day */
    kv = exp(7.755 + 0.0352 * silt + 0.93 * (double)topsoil -
        0.967 * bd * bd - 0.000484 * clay * clay - 0.000322 * silt * silt +
        0.001 / silt - 0.0748 / om - 0.643 * log(silt) -
        0.01398 * bd * clay - 0.1673 * bd * om +
        0.02986 * (double)topsoil * clay - 0.03305 * (double)topsoil * silt);

    /* Convert from cm/day to m/s */
    kv /= 100.0 * 24.0 * 3600.0;

    return kv;
}

double PtfThetas(double silt, double clay, double om, double bd, int topsoil)
{
    double          thetas;

    thetas = 0.7919 + 0.001691 * clay - 0.29619 * bd -
        0.000001491 * silt * silt + 0.0000821 * om * om + 0.02427 / clay +
        0.01113 / silt + 0.01472 * log(silt) - 0.0000733 * om * clay -
        0.000619 * bd * clay - 0.001183 * bd * om -
        0.0001664 * (double)topsoil *silt;

    return thetas;
}

double PtfThetar(double silt, double clay)
{
    double          thetar;

    if (clay < 18.0 && 100.0 - silt - clay > 65.0)
    {
        thetar = 0.15;
    }
    else
    {
        thetar = 0.05;
    }

    return thetar;
}

double PtfAlpha(double silt, double clay, double om, double bd, int topsoil)
{
    double          alpha;

    /* Calcualte alpha in cm */
    alpha = exp(-14.96 + 0.03135 * clay + 0.0351 * silt + 0.646 * om +
        15.29 * bd - 0.192 * (double)topsoil - 4.671 * bd * bd -
        0.000781 * clay * clay - 0.00687 * om * om + 0.0449 / om +
        0.0663 * log(silt) + 0.1482 * log(om) - 0.04546 * bd * silt -
        0.4852 * bd * om + 0.00673 * (double)topsoil * clay);

    /* Convert from 1/cm to 1/m */
    alpha *= 100.0;

    return alpha;
}

double PtfBeta(double silt, double clay, double om, double bd, int topsoil)
{
    double          beta;

    beta = 1.0 + exp(-25.23 - 0.02195 * clay + 0.0074 * silt - 0.1940 * om +
        45.5 * bd - 7.24 * bd * bd + 0.0003658 * clay * clay +
        0.002885 * om * om - 12.81 / bd - 0.1524 / silt - 0.01958 / om -
        0.2876 * log(silt) - 0.0709 * log(om) - 44.6 * log(bd) -
        0.02264 * bd * clay + 0.0896 * bd * om +
        0.00718 * (double)topsoil * clay);

    return beta;
}
