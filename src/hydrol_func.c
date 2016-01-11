#include "pihm.h"

double DhByDl (double *l1, double *l2, double *surfh)
{
    return (-1.0 * (l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] -
                surfh[2]) + l1[0] * (surfh[2] - surfh[1])) / (l2[2] * (l1[1] -
                l1[0]) + l2[1] * (l1[0] - l1[2]) + l2[0] * (l1[2] - l1[1])));
}

double RivArea (int riv_order, double riv_depth, double riv_coeff)
{
    double          riv_area;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            riv_area = riv_depth * riv_coeff;
            break;
        case TRIANGLE:
            riv_area = pow (riv_depth, 2) / riv_coeff;
            break;
        case QUADRATIC:
            riv_area =
                4.0 * pow (riv_depth, 1.5) / (3.0 * pow (riv_coeff, 0.5));
            break;
        case CUBIC:
            riv_area =
                3.0 * pow (riv_depth, 4.0 / 3.0) / (2.0 * pow (riv_coeff,
                    1.0 / 3.0));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }

    return (riv_area);
}

double RivPerim (int riv_order, double riv_depth, double riv_coeff)
{
    double          riv_perim;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            riv_perim = 2.0 * riv_depth + riv_coeff;
            break;
        case TRIANGLE:
            riv_perim =
                2.0 * riv_depth * pow (1.0 + pow (riv_coeff, 2),
                0.5) / riv_coeff;
            break;
        case QUADRATIC:
            riv_perim =
                (pow (riv_depth * (1.0 +
                        4.0 * riv_coeff * riv_depth) / riv_coeff,
                    0.5)) + (log (2.0 * pow (riv_coeff * riv_depth,
                        0.5) + pow (1.0 + 4.0 * riv_coeff * riv_depth,
                        0.5)) / (2.0 * riv_coeff));
            break;
        case CUBIC:
            riv_perim =
                2.0 * ((pow (riv_depth * (1.0 + 9.0 * pow (riv_coeff,
                                2.0 / 3.0) * riv_depth),
                        0.5) / 3.0) + (log (3.0 * pow (riv_coeff,
                            1.0 / 3.0) * pow (riv_depth,
                            0.5) + pow (1.0 + 9.0 * pow (riv_coeff,
                                2.0 / 3.0) * riv_depth,
                            0.5)) / (9.0 * pow (riv_coeff, 1.0 / 3.0))));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }
    return (riv_perim);
}

double EqWid (int riv_order, double riv_depth, double riv_coeff)
{
    double          eq_wid;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            eq_wid = riv_coeff;
            break;
        case TRIANGLE:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        case QUADRATIC:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        case CUBIC:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }
    return (eq_wid);
}


double OLFEleToRiv (double eleytot, double elez, double cwr, double rivzmax,
    double rivytot, double length)
{
    double          flux;
    double          zbank;

    /*
     * Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
     */
    zbank = (rivzmax < elez) ? elez : rivzmax;

    if (rivytot > eleytot)
    {
        if (eleytot > zbank)
        {
            /* Submerged weir */
            flux =
                cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (rivytot -
                eleytot) * (rivytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < rivytot)
            {
                /* Free-flowing weir */
                flux =
                    cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (rivytot -
                    zbank) * (rivytot - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else
    {
        if (rivytot > zbank)
        {
            /* Submerged weir */
            flux =
                -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleytot -
                rivytot) * (eleytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < eleytot)
            {
                /* Free-flowing weir */
                flux =
                    -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleytot -
                    zbank) * (eleytot - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }

    return (flux);
}

double OverlandFlow (double avg_y, double grad_y, double avg_sf,
    double crossa, double avg_rough)
{
    return (crossa * pow (avg_y,
            2.0 / 3.0) * grad_y / (sqrt (fabs (avg_sf)) * avg_rough));
}

double AvgY (double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > 1.0 * IMMOBILE)
        {
            avg_y = 1.0 * yi;
        }
        else
        {
            avg_y = 0.0;
        }
    }
    else
    {
        if (yinabr > 1.0 * IMMOBILE)
        {
            avg_y = 1.0 * yinabr;
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return (avg_y);
}

double EffKV (double ksatfunc, double elemsatn, int status, double mackv,
    double kv, double areaf)
{
    double          effkv;

    switch (status)
    {
        case MAT_CTRL:
            effkv = kv * ksatfunc;
            break;
        case APP_CTRL:
            effkv = mackv * areaf * elemsatn + kv * (1.0 - areaf) * ksatfunc;
            break;
        case MAC_CTRL:
            effkv = mackv * areaf + kv * (1.0 - areaf) * ksatfunc;
            break;
        case SAT_CTRL:
            effkv = mackv * areaf + kv * (1.0 - areaf) * ksatfunc;
            break;
        default:
            printf ("Error: Macropore status not recognized!\n");
            PihmExit (1);
    }

    return (effkv);
}

int MacroporeStatus (double ksatfunc, double elemsatn, double grady,
    double mackv, double kv, double areaf)
{
    if (elemsatn >= 0.98)
    {
        return (SAT_CTRL);
    }
    else
    {
        if (grady * ksatfunc * kv <= 1.0 * kv * ksatfunc)
        {
            return (MAT_CTRL);
        }
        else
        {
            if (grady * ksatfunc * kv <
                (mackv * areaf + kv * (1.0 - areaf) * ksatfunc))
            {
                return (APP_CTRL);
            }
            else
            {
                return (MAC_CTRL);
            }
        }
    }
}

double EffKH (int mp, double tmpy, double aqdepth, double macd,
    double macksath, double areaf, double ksath)
{
    tmpy = (tmpy > 0.0) ? tmpy : 0.0;

    if (mp == 1)
    {
        if (tmpy > aqdepth - macd)
        {
            if (tmpy > aqdepth)
            {
                return (macksath * macd * areaf + ksath * (aqdepth -
                        macd * areaf)) / aqdepth;
            }
            else
            {
                return (macksath * (tmpy - (aqdepth - macd)) * areaf +
                    ksath * (aqdepth - macd + (tmpy - (aqdepth -
                                macd)) * (1.0 - areaf))) / tmpy;
            }
        }
        else
        {
            return (ksath);
        }
    }
    else
    {
        return (ksath);
    }
}

double KrFunc (double alpha, double beta, double satn)
{
    return (pow (satn, 0.5) * pow (-1.0 + pow (1.0 - pow (satn,
                    beta / (beta - 1.0)), (beta - 1.0) / beta), 2));
}

double Psi (double satn, double alpha, double beta)
{
    /* van Genuchten 1980 SSSAJ */
    return (0.0 - pow (pow (1.0 / satn, beta / (beta - 1.0)) - 1.0,
            1.0 / beta) / alpha);
}
