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
        elem->wf.prcp = *elem->forc.meteo[PRCP_TS] / 1000.0;
        elem->ps.albedo = 0.5 * (elem->lc.albedomin + elem->lc.albedomax);
        radnet = *elem->forc.meteo[SOLAR_TS] * (1.0 - elem->ps.albedo);
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
        snow_rate = frac_snow * elem->wf.prcp;
        /* snow_grnd, snow_canopy, snow_intcp_max, melt_rate_grnd,
         * melt_rate_canopy are the average value prorated over the whole
         * elemental area */
        snow_grnd =
            snow_grnd + (1.0 - elem->lc.shdfac) * snow_rate * stepsize;
        snow_canopy = snow_canopy + elem->lc.shdfac * snow_rate * stepsize;
        snow_intcp_max =
            snow_canopy > 0.0 ? 0.003 * lai * elem->lc.shdfac : 0.0;
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
        intcp_max = elem->lc.cmcfactr * lai * elem->lc.shdfac;

        z0 = MonthlyRL (t, elem->lc.type);

        ra = log (elem->ps.zlvl_wind / z0) * log (10.0 *
            elem->ps.zlvl_wind / z0) / (wind * 0.16);

        gamma =
            4.0 * 0.7 * SIGMA * RD / CP * pow (sfctmp + 273.15,
            4) / (pres / ra) + 1.0;
        delta =
            LVH2O * LVH2O * 0.622 / RV / CP / pow (sfctmp + 273.15,
            2) * qvsat;

        etp =
            (radnet * delta + gamma * (1.2 * LVH2O * (qvsat -
                    qv) / ra)) / (1000.0 * LVH2O * (delta + gamma));

        if (elem->soil.depth - elem->ws0.gw < elem->lc.rzd)
            satn = 1.0;
        else
            satn =
                ((elem->ws0.unsat / (elem->soil.depth - elem->ws0.gw)) >
                1.0) ? 1.0 : ((elem->ws0.unsat / (elem->soil.depth -
                        elem->ws0.gw)) <
                0.0) ? 0.0 : 0.5 * (1.0 -
                cos (3.14 * (elem->ws0.unsat / (elem->soil.depth - elem->ws0.gw))));

        betas =
            (satn * elem->soil.porosity + elem->soil.smcmin -
            elem->soil.smcwlt) / (elem->soil.smcref - elem->soil.smcwlt);
        betas = (betas < 0.0001) ? 0.0001 : (betas > 1.0 ? 1.0 : betas);
        elem->wf.edir = (1.0 - elem->lc.shdfac) * pow (betas, 2) * etp;
        elem->wf.edir = elem->wf.edir < 0.0 ? 0.0 : elem->wf.edir;

        /* Note the dependence on physical units */
        if (lai > 0.0)
        {
            elem->wf.ec =
                elem->lc.shdfac * (pow ((elem->ws.cmc <
                        0.0 ? 0.0 : (elem->ws.cmc >
                            intcp_max ? intcp_max : elem->ws.cmc)) / intcp_max,
                    elem->lc.cfactr)) * etp;
            elem->wf.ec = elem->wf.ec < 0.0 ? 0.0 : elem->wf.ec;

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

            elem->wf.ett =
                elem->lc.shdfac * pc * (1.0 - pow (((elem->ws.cmc +
                            snow_canopy <
                            0.0) ? 0.0 : (elem->ws.cmc +
                            snow_canopy)) / (intcp_max + snow_intcp_max),
                    elem->lc.cfactr)) * etp;
            elem->wf.ett = elem->wf.ett < 0.0 ? 0.0 : elem->wf.ett;
            elem->wf.ett = ((elem->ws.gw < (elem->soil.depth - elem->lc.rzd)) &&
                elem->ws.unsat <= 0.0) ? 0.0 : elem->wf.ett;

            elem->wf.drip =
                elem->ws.cmc <=
                0.0 ? 0.0 : 5.65e-2 * intcp_max * exp (3.89 * (elem->ws.cmc <
                    0.0 ? 0.0 : elem->ws.cmc) / intcp_max);
        }
        else
        {
            elem->wf.ett = 0.0;
            elem->wf.ec = 0.0;
            elem->wf.drip = 0.0;
        }

        if (elem->wf.drip < 0.0)
            elem->wf.drip = 0.0;
        if (elem->wf.drip * stepsize > elem->ws.cmc)
            elem->wf.drip = elem->ws.cmc / stepsize;
        if (elem->ws.cmc >= intcp_max)
        {
            if (((1.0 - frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                    melt_rate_canopy) >= elem->wf.ec + elem->wf.drip)
            {
                ret =
                    elem->wf.drip + (elem->ws.cmc - intcp_max) / stepsize +
                    (((1.0 - frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy) - (elem->wf.ec + elem->wf.drip));
                isval = intcp_max;
            }
            else if ((((1.0 - frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy) < elem->wf.ec + elem->wf.drip) &&
                (elem->ws.cmc + stepsize * ((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy - elem->wf.ec - elem->wf.drip) <= 0.0))
            {
                elem->wf.ec =
                    (elem->wf.ec / (elem->wf.ec +
                        elem->wf.drip)) * (elem->ws.cmc / stepsize + ((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy));
                ret =
                    (elem->wf.drip / (elem->wf.ec +
                        elem->wf.drip)) * (elem->ws.cmc / stepsize + ((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy));
                isval = 0.0;
            }
            else
            {
                isval =
                    elem->ws.cmc + stepsize * (((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy) - elem->wf.ec - elem->wf.drip);
                ret = elem->wf.drip;
            }
        }
        else if ((elem->ws.cmc < intcp_max) &&
            ((elem->ws.cmc + (((1.0 -
                                frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                            melt_rate_canopy) - elem->wf.ec -
                        elem->wf.drip) * stepsize) >= intcp_max))
        {
            isval = intcp_max;
            ret =
                elem->wf.drip + (((elem->ws.cmc + (((1.0 -
                                    frac_snow) * elem->wf.prcp *
                                elem->lc.shdfac + melt_rate_canopy) -
                            elem->wf.ec - elem->wf.drip) * stepsize) -
                    intcp_max)) / stepsize;
        }
        else if ((elem->ws.cmc < intcp_max) &&
            ((elem->ws.cmc + (((1.0 -
                                frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                            melt_rate_canopy) - elem->wf.ec -
                        elem->wf.drip) * stepsize) <= 0.0))
        {
            if ((elem->wf.ec > 0.0) || (elem->wf.drip > 0.0))
            {
                elem->wf.ec =
                    (elem->wf.ec / (elem->wf.ec +
                        elem->wf.drip)) * (elem->ws.cmc / stepsize + ((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy));
                ret =
                    (elem->wf.drip / (elem->wf.ec +
                        elem->wf.drip)) * (elem->ws.cmc / stepsize + ((1.0 -
                            frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                        melt_rate_canopy));
            }
            else
            {
                elem->wf.ec = 0.0;
                ret = 0.0;
            }
            isval = 0.0;
        }
        else
        {
            isval =
                elem->ws.cmc + (((1.0 -
                        frac_snow) * elem->wf.prcp * elem->lc.shdfac +
                    melt_rate_canopy) - elem->wf.ec - elem->wf.drip) * stepsize;
            ret = elem->wf.drip;
        }

        elem->wf.netprcp =
            (1.0 - elem->lc.shdfac) * (1.0 - frac_snow) * elem->wf.prcp + ret +
            melt_rate_grnd;
        elem->wf.drip = ret;
        elem->ws.cmc = isval;
    }
}
