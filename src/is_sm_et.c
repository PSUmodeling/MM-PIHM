#include "pihm.h"

#if !defined(_NOAH_)
void IntcpSnowEt(int t, double stepsize, elem_struct *elem,
    const calib_struct *cal)
{
    int             i;
    const double    TSNOW = -3.0;
    const double    TRAIN = 1.0;
    const double    T0 = 0.0;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
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
        double          melt_rate;
        double          intcp_max;
        double          meltf;

        /* Note the dependence on physical units */
        elem[i].ps.albedo = 0.5 * (elem[i].lc.albedomin + elem[i].lc.albedomax);
        radnet = elem[i].ef.soldn * (1.0 - elem[i].ps.albedo);
        sfctmp = elem[i].es.sfctmp - 273.15;
        wind = elem[i].ps.sfcspd;
        rh = elem[i].ps.rh / 100.0;

        vp = 611.2 * exp(17.67 * sfctmp / (sfctmp + 243.5)) * rh;
        pres = 101.325 * 1.0e3 *
            pow((293.0 - 0.0065 * elem[i].topo.zmax) / 293.0, 5.26);
        qv = 0.622 * vp / pres;
        qvsat = 0.622 * (vp / rh) / pres;

        if (elem[i].attrib.lai_type > 0)
        {
            lai = elem[i].ps.proj_lai;
        }
        else
        {
            lai = MonthlyLai(t, elem[i].attrib.lc_type);
        }

        meltf = MonthlyMf(t);

        /* Snow accumulation and snow melt calculation */
        frac_snow = (sfctmp < TSNOW) ?
            1.0 :
            ((sfctmp > TRAIN) ? 0.0 : (TRAIN - sfctmp) / (TRAIN - TSNOW));
        snow_rate = frac_snow * elem[i].wf.prcp;

        elem[i].ws.sneqv += snow_rate * stepsize;

        melt_rate = (sfctmp > T0) ? (sfctmp - T0) * meltf : 0.0;

        if (elem[i].ws.sneqv > melt_rate * stepsize)
        {
            elem[i].ws.sneqv -= melt_rate * stepsize;
        }
        else
        {
            melt_rate = elem[i].ws.sneqv / stepsize;
            elem[i].ws.sneqv = 0.0;
        }

        /* ThroughFall and Evaporation from canopy */
        intcp_max = elem[i].lc.cmcfactr * lai * elem[i].lc.shdfac;

        z0 = MonthlyRl(t, elem[i].attrib.lc_type);

        ra = log(elem[i].ps.zlvl_wind / z0) *
            log(10.0 * elem[i].ps.zlvl_wind / z0) / (wind * 0.16);

        gamma = 4.0 * 0.7 * SIGMA * RD / CP * pow(sfctmp + 273.15, 4) /
            (pres / ra) + 1.0;
        delta =
            LVH2O * LVH2O * 0.622 / RV / CP / pow(sfctmp + 273.15, 2) * qvsat;

        etp = (radnet * delta + gamma * (1.2 * LVH2O * (qvsat - qv) / ra)) /
            (1000.0 * LVH2O * (delta + gamma));

        if (elem[i].soil.depth - elem[i].ws.gw < elem[i].ps.rzd)
        {
            satn = 1.0;
        }
        else
        {
            satn =
                ((elem[i].ws.unsat /
                (elem[i].soil.depth - elem[i].ws.gw)) > 1.0) ?
                1.0 :
                ((elem[i].ws.unsat /
                (elem[i].soil.depth - elem[i].ws.gw)) < 0.0) ?
                0.0 :
                0.5 * (1.0 - cos(3.14 *
                (elem[i].ws.unsat / (elem[i].soil.depth - elem[i].ws.gw))));
        }

        betas = (satn * elem[i].soil.porosity +
            elem[i].soil.smcmin - elem[i].soil.smcwlt) /
            (elem[i].soil.smcref - elem[i].soil.smcwlt);
        betas = (betas < 0.0001) ? 0.0001 : ((betas > 1.0) ? 1.0 : betas);
        elem[i].wf.edir = (1.0 - elem[i].lc.shdfac) * pow(betas, 2) * etp;
        elem[i].wf.edir *= cal->edir;
        elem[i].wf.edir = (elem[i].wf.edir < 0.0) ? 0.0 : elem[i].wf.edir;

        /* Note the dependence on physical units */
        if (lai > 0.0)
        {
            elem[i].wf.ec = elem[i].lc.shdfac *
                pow(((elem[i].ws.cmc < 0.0) ? 0.0 :
                ((elem[i].ws.cmc > intcp_max) ? intcp_max : elem[i].ws.cmc)) /
                intcp_max, elem[i].lc.cfactr) * etp;
            elem[i].wf.ec *= cal->ec;
            elem[i].wf.ec = (elem[i].wf.ec < 0.0) ? 0.0 : elem[i].wf.ec;

            fr = 1.1 * radnet / (elem[i].epc.rgl * lai);
            fr = (fr < 0.0) ? 0.0 : fr;
            alphar =
                (1.0 + fr) / (fr + (elem[i].epc.rsmin / elem[i].epc.rsmax));
            alphar = (alphar > 10000.0) ? 10000.0 : alphar;
            etas =
                1.0 - 0.0016 * (pow((elem[i].epc.topt - 273.15 - sfctmp), 2));
            etas = (etas < 0.0001) ? 0.0001 : etas;
            gammas = 1.0 / (1.0 + 0.00025 * (vp / rh - vp));
            gammas = (gammas < 0.01) ? 0.01 : gammas;
            rs = elem[i].epc.rsmin * alphar / (betas * lai * etas * gammas);
            rs = (rs > elem[i].epc.rsmax) ? elem[i].epc.rsmax : rs;

            pc = (1.0 + delta / gamma) / (1.0 + rs / ra + delta / gamma);

            elem[i].wf.ett = elem[i].lc.shdfac * pc *
                (1.0 - pow((elem[i].ws.cmc < 0.0) ?
                0.0 :
                ((elem[i].ws.cmc > intcp_max) ?
                intcp_max : elem[i].ws.cmc) / intcp_max, elem[i].lc.cfactr)) *
                etp;
            elem[i].wf.ett *= cal->ett;
            elem[i].wf.ett = (elem[i].wf.ett < 0.0) ? 0.0 : elem[i].wf.ett;
            elem[i].wf.ett =
                ((elem[i].ws.gw < (elem[i].soil.depth - elem[i].ps.rzd))
                && elem[i].ws.unsat <= 0.0) ? 0.0 : elem[i].wf.ett;

            /* Drip function from Rutter and Morton, 1977, Journal of Applied
             * Ecology
             * D0 = 3.91E-5 m/min = 6.52E-7 m/s */
            elem[i].wf.drip = (elem[i].ws.cmc <= 0.0) ?
                0.0 :
                6.52E-7 * intcp_max * exp(3.89 * elem[i].ws.cmc / intcp_max);
        }
        else
        {
            elem[i].wf.ett = 0.0;
            elem[i].wf.ec = 0.0;
            elem[i].wf.drip = 0.0;
        }

        if (elem[i].wf.drip < 0.0)
        {
            elem[i].wf.drip = 0.0;
        }
        if (elem[i].wf.drip * stepsize > elem[i].ws.cmc)
        {
            elem[i].wf.drip = elem[i].ws.cmc / stepsize;
        }

        isval = elem[i].ws.cmc +
            (1.0 - frac_snow) * elem[i].wf.prcp * elem[i].lc.shdfac * stepsize -
            elem[i].wf.ec * stepsize - elem[i].wf.drip * stepsize;

        if (isval > intcp_max)
        {
            elem[i].ws.cmc = intcp_max;
            elem[i].wf.drip += (isval - intcp_max) / stepsize;
        }
        else if (isval < 0.0)
        {
            elem[i].ws.cmc = 0.0;
            if (elem[i].wf.ec + elem[i].wf.drip > 0.0)
            {
                elem[i].wf.ec =
                    elem[i].wf.ec / (elem[i].wf.ec + elem[i].wf.drip) *
                    (elem[i].ws.cmc +
                    (1.0 - frac_snow) * elem[i].wf.prcp * elem[i].lc.shdfac *
                    stepsize);
                elem[i].wf.drip =
                    elem[i].wf.drip / (elem[i].wf.ec + elem[i].wf.drip) *
                    (elem[i].ws.cmc +
                    (1.0 - frac_snow) * elem[i].wf.prcp * elem[i].lc.shdfac *
                    stepsize);
            }
        }
        else
        {
            elem[i].ws.cmc = isval;
        }

        elem[i].wf.pcpdrp =
            (1.0 - elem[i].lc.shdfac) * (1.0 - frac_snow) * elem[i].wf.prcp +
            elem[i].wf.drip + melt_rate;
    }
}
#endif
