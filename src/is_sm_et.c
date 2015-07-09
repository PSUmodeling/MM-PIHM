#include "pihm.h"

void IntcpSnowET (int t, double stepsize, pihm_struct pihm)
{
    int             i;
    double          satn;
    double          betas;
    double          fr;
    double          alphar;
    double          etas;
    double          gammas;
    double          rs;
    double          pc;
    double          delta, gamma;
    double          radnet;
    double          sfctmp;
    double          wind;
    double          rh;
    double          vp;
    double          pres;
    double          lai;
    double          z0;
    double          ra;
    double          qvsat;
    double          qv;
    double          etp;
    double          isval = 0.0;
    double          frac_snow;
    double          snow_rate;
    double          melt_rate_grnd;
    double          melt_rate_canopy;
    double          snow_grnd = 0.0;
    double          snow_canopy = 0.0;
    double          snow_intcp_max;
    double          intcp_max;
    double          meltf;
    const double    TSNOW = -3.0;
    const double    TRAIN = 1.0;
    const double    T0 = 0.0;
    double          ret = 0.0;

    elem_struct    *elem;

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        /* Note the dependence on physical units */
        elem->prcp = *elem->forc.meteo[PRCP_TS] / 1000.0;
        elem->albedo = 0.5 * (elem->lc.albedomin + elem->lc.albedomax);
        radnet = *elem->forc.meteo[SOLAR_TS] * (1.0 - elem->albedo);
        sfctmp = *elem->forc.meteo[SFCTMP_TS] - 273.15;
        wind = *elem->forc.meteo[SFCSPD_TS];
        rh = *elem->forc.meteo[RH_TS] / 100.0;

        vp = 611.2 * exp (17.67 * sfctmp / (sfctmp + 243.5)) * rh;
        pres =
            101.325 * 1.0e3 * pow ((293.0 - 0.0065 * elem->topo.zmax) / 293.0,
            5.26);
        qv = 0.622 * vp / pres;
        qvsat = 0.622 * (vp / rh) / pres;
        if (elem->forc.lai_type > 0)
        {
            lai = *elem->forc.lai;
        }
        else
        {
            lai = MonthlyLAI (t, elem->lc.type);
        }

        meltf = MonthlyMF (t);

        /* Snow accumulation and snow melt calculation */
        frac_snow =
            sfctmp < TSNOW ? 1.0 : sfctmp >
            TRAIN ? 0 : (TRAIN - sfctmp) / (TRAIN - TSNOW);
        snow_rate = frac_snow * elem->prcp;
        /* snow_grnd, snow_canopy, snow_intcp_max, melt_rate_grnd,
         * melt_rate_canopy are the average value prorated over the whole
         * elemental area */
        snow_grnd =
            snow_grnd + (1.0 - elem->lc.vegfrac) * snow_rate * stepsize;
        snow_canopy = snow_canopy + elem->lc.vegfrac * snow_rate * stepsize;
        snow_intcp_max =
            snow_canopy > 0.0 ? 0.003 * lai * elem->lc.vegfrac : 0.0;
        if (snow_canopy > snow_intcp_max)
        {
            snow_grnd += snow_canopy - snow_intcp_max;
            snow_canopy = snow_intcp_max;
        }
        melt_rate_grnd = melt_rate_canopy =
            (sfctmp > T0 ? (sfctmp - T0) * meltf : 0.0);

        if (snow_grnd > melt_rate_grnd * stepsize)
            snow_grnd = snow_grnd - melt_rate_grnd * stepsize;
        else
        {
            melt_rate_grnd = snow_grnd / stepsize;
            snow_grnd = 0.0;
        }
        if (snow_canopy > melt_rate_canopy * stepsize)
            snow_canopy = snow_canopy - melt_rate_canopy * stepsize;
        else
        {
            melt_rate_canopy = snow_canopy / stepsize;
            snow_canopy = 0.0;
        }

        /* ThroughFall and Evaporation from canopy */
        /*
         * EleIS, EleET[0] and ret are prorated for the whole element.
         * Logistics are simpler if assumed in volumetric form by
         * multiplication of Area on either side of equation */
        intcp_max = elem->lc.intcp_factr * lai * elem->lc.vegfrac;

        z0 = MonthlyRL (t, elem->lc.type);

        ra = log (elem->forc.zlvl_wind / z0) * log (10.0 *
            elem->forc.zlvl_wind / z0) / (wind * 0.16);

        gamma =
            4.0 * 0.7 * SIGMA * RD / CP * pow (sfctmp + 273.15,
            4) / (pres / ra) + 1.0;
        delta =
            LVH2O * LVH2O * 0.622 / RV / CP / pow (sfctmp + 273.15,
            2) * qvsat;

        etp =
            (radnet * delta + gamma * (1.2 * LVH2O * (qvsat -
                    qv) / ra)) / (1000.0 * LVH2O * (delta + gamma));

        if (elem->soil.depth - elem->gw0 < elem->lc.rzd)
            satn = 1.0;
        else
            satn =
                ((elem->unsat0 / (elem->soil.depth - elem->gw0)) >
                1.0) ? 1.0 : ((elem->unsat0 / (elem->soil.depth -
                        elem->gw0)) <
                0.0) ? 0.0 : 0.5 * (1.0 -
                cos (3.14 * (elem->unsat0 / (elem->soil.depth - elem->gw0))));

        betas =
            (satn * elem->soil.porosity + elem->soil.thetar -
            elem->soil.thetaw) / (elem->soil.thetaref - elem->soil.thetaw);
        betas = (betas < 0.0001) ? 0.0001 : (betas > 1.0 ? 1.0 : betas);
        elem->et[2] = (1.0 - elem->lc.vegfrac) * pow (betas, 2) * etp;
        elem->et[2] = elem->et[2] < 0.0 ? 0.0 : elem->et[2];

        /* Note the dependence on physical units */
        if (lai > 0.0)
        {
            elem->et[0] =
                elem->lc.vegfrac * (pow ((elem->intcp <
                        0.0 ? 0.0 : (elem->intcp >
                            intcp_max ? intcp_max : elem->intcp)) / intcp_max,
                    elem->lc.cfactr)) * etp;
            elem->et[0] = elem->et[0] < 0.0 ? 0.0 : elem->et[0];

            fr = 1.1 * radnet / (elem->lc.rgl * lai);
            fr = fr < 0.0 ? 0.0 : fr;
            alphar = (1.0 + fr) / (fr + (elem->lc.rsmin / elem->lc.rsmax));
            alphar = alphar > 10000.0 ? 10000.0 : alphar;
            etas =
                1.0 - 0.0016 * (pow ((elem->lc.topt - 273.15 - sfctmp), 2));
            etas = etas < 0.0001 ? 0.0001 : etas;
            gammas = 1.0 / (1.0 + 0.00025 * (vp / rh - vp));
            gammas = (gammas < 0.01) ? 0.01 : gammas;
            rs = elem->lc.rsmin * alphar / (betas * lai * etas * gammas);
            rs = rs > elem->lc.rsmax ? elem->lc.rsmax : rs;

            pc = (1.0 + delta / gamma) / (1.0 + rs / ra + delta / gamma);

            elem->et[1] =
                elem->lc.vegfrac * pc * (1.0 - pow (((elem->intcp +
                            snow_canopy <
                            0.0) ? 0.0 : (elem->intcp +
                            snow_canopy)) / (intcp_max + snow_intcp_max),
                    elem->lc.cfactr)) * etp;
            elem->et[1] = elem->et[1] < 0.0 ? 0.0 : elem->et[1];
            elem->et[1] = ((elem->gw < (elem->soil.depth - elem->lc.rzd)) &&
                elem->unsat <= 0.0) ? 0.0 : elem->et[1];

            elem->drip =
                elem->intcp <=
                0.0 ? 0.0 : 5.65e-2 * intcp_max * exp (3.89 * (elem->intcp <
                    0.0 ? 0.0 : elem->intcp) / intcp_max);
        }
        else
        {
            elem->et[1] = 0.0;
            elem->et[0] = 0.0;
            elem->drip = 0.0;
        }

        if (elem->drip < 0.0)
            elem->drip = 0.0;
        if (elem->drip * stepsize > elem->intcp)
            elem->drip = elem->intcp / stepsize;
        if (elem->intcp >= intcp_max)
        {
            if (((1.0 - frac_snow) * elem->prcp * elem->lc.vegfrac +
                    melt_rate_canopy) >= elem->et[0] + elem->drip)
            {
                ret =
                    elem->drip + (elem->intcp - intcp_max) / stepsize +
                    (((1.0 - frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy) - (elem->et[0] + elem->drip));
                isval = intcp_max;
            }
            else if ((((1.0 - frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy) < elem->et[0] + elem->drip) &&
                (elem->intcp + stepsize * ((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy - elem->et[0] - elem->drip) <= 0.0))
            {
                elem->et[0] =
                    (elem->et[0] / (elem->et[0] +
                        elem->drip)) * (elem->intcp / stepsize + ((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy));
                ret =
                    (elem->drip / (elem->et[0] +
                        elem->drip)) * (elem->intcp / stepsize + ((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy));
                isval = 0.0;
            }
            else
            {
                isval =
                    elem->intcp + stepsize * (((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy) - elem->et[0] - elem->drip);
                ret = elem->drip;
            }
        }
        else if ((elem->intcp < intcp_max) &&
            ((elem->intcp + (((1.0 -
                                frac_snow) * elem->prcp * elem->lc.vegfrac +
                            melt_rate_canopy) - elem->et[0] -
                        elem->drip) * stepsize) >= intcp_max))
        {
            isval = intcp_max;
            ret =
                elem->drip + (((elem->intcp + (((1.0 -
                                    frac_snow) * elem->prcp *
                                elem->lc.vegfrac + melt_rate_canopy) -
                            elem->et[0] - elem->drip) * stepsize) -
                    intcp_max)) / stepsize;
        }
        else if ((elem->intcp < intcp_max) &&
            ((elem->intcp + (((1.0 -
                                frac_snow) * elem->prcp * elem->lc.vegfrac +
                            melt_rate_canopy) - elem->et[0] -
                        elem->drip) * stepsize) <= 0.0))
        {
            if ((elem->et[0] > 0.0) || (elem->drip > 0.0))
            {
                elem->et[0] =
                    (elem->et[0] / (elem->et[0] +
                        elem->drip)) * (elem->intcp / stepsize + ((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy));
                ret =
                    (elem->drip / (elem->et[0] +
                        elem->drip)) * (elem->intcp / stepsize + ((1.0 -
                            frac_snow) * elem->prcp * elem->lc.vegfrac +
                        melt_rate_canopy));
            }
            else
            {
                elem->et[0] = 0.0;
                ret = 0.0;
            }
            isval = 0.0;
        }
        else
        {
            isval =
                elem->intcp + (((1.0 -
                        frac_snow) * elem->prcp * elem->lc.vegfrac +
                    melt_rate_canopy) - elem->et[0] - elem->drip) * stepsize;
            ret = elem->drip;
        }

        elem->netprcp =
            (1.0 - elem->lc.vegfrac) * (1.0 - frac_snow) * elem->prcp + ret +
            melt_rate_grnd;
        elem->drip = ret;
        elem->intcp = isval;
    }
}
