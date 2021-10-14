#include "pihm.h"

void VerticalFlow(double dt, elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        // Calculate infiltration rate
        elem[i].wf.infil = Infil(dt, &elem[i].topo, &elem[i].soil, &elem[i].ws, &elem[i].ws0, &elem[i].wf);

#if defined(_NOAH_)
        // Constrain infiltration by frozen top soil
        elem[i].wf.infil *= elem[i].ps.fcr;
#endif

        // Calculate recharge rate
        elem[i].wf.recharge = Recharge(&elem[i].soil, &elem[i].ws, &elem[i].wf);

#if defined(_DGW_)
        elem[i].wf.infil_geol = GeolInfil(&elem[i].topo, &elem[i].soil, &elem[i].geol, &elem[i].ws);
        elem[i].wf.rechg_geol = GeolRecharge(&elem[i].geol, &elem[i].ws, &elem[i].wf);
#endif
    }
}

double Infil(double dt, const topo_struct *topo, const soil_struct *soil, const wstate_struct *ws,
    const wstate_struct *ws0, const wflux_struct *wf)
{
    double          appl_rate;
    double          wet_frac;
    double          dh_dz;
    double          satn;
    double          satkfunc;
    double          infil;
    double          infil_max;
    double          kinf;
    double          deficit;
    double          psi_u;
    double          h_u;
    int             j;

    if (ws->gw > soil->depth)
    {
        infil = -soil->kinfv;
    }
    else if (ws->unsat + ws->gw > soil->depth)
    {
        infil = 0.0;
    }
    else
    {
        // Water application rate is the sum of net gain of overland flow and net precipitation
        appl_rate = 0.0;
        for (j = 0; j < NUM_EDGE; j++)
        {
            appl_rate += -wf->overland[j] / topo->area;
        }
        appl_rate = MAX(appl_rate, 0.0);
        appl_rate += wf->pcpdrp;

        if (DEPRSTG > 0.0)
        {
            wet_frac = ws->surfh / DEPRSTG;
            wet_frac = MAX(wet_frac, 0.0);
            wet_frac = MIN(wet_frac, 1.0);
        }
        else
        {
            wet_frac = (ws->surfh > 0.0) ? 1.0 : 0.0;
        }

        if (ws->gw > soil->depth - soil->dinf)
        {
            // Assumption: Dinf < Dmac
            dh_dz = (ws->surfh + topo->zmax - (ws->gw + topo->zmin)) / (0.5 * (ws->surfh + soil->dinf));
            dh_dz = (ws->surfh <= 0.0 && dh_dz > 0.0) ? 0.0 : dh_dz;

            satn = 1.0;
            satkfunc = KrFunc(soil->beta, satn);

            kinf = EffKinf(dh_dz, satkfunc, satn, appl_rate, ws->surfh, soil);

            infil = kinf * dh_dz;
        }
        else
        {
            deficit = soil->depth - ws->gw;
#if defined(_NOAH_)
            satn = (ws->swc[0] - soil->smcmin) / (soil->smcmax - soil->smcmin);
#else
            satn = ws->unsat / deficit;
#endif
            satn = MIN(satn, 1.0);
            satn = MAX(satn, SATMIN);

            psi_u = Psi(satn, soil->alpha, soil->beta);
            // Note: for psi calculation using van Genuchten relation, cutting the psi-sat tail at small saturation can
            // be performed for computational advantage. If you do not want to perform this, comment the statement that
            // follows
            psi_u = MAX(psi_u, PSIMIN);

            h_u = psi_u + topo->zmax - 0.5 * soil->dinf;
            dh_dz = (ws->surfh + topo->zmax - h_u) / (0.5 * (ws->surfh + soil->dinf));
            dh_dz = (ws->surfh <= 0.0 && dh_dz > 0.0) ?  0.0 : dh_dz;

            satkfunc = KrFunc(soil->beta, satn);

            kinf = EffKinf(dh_dz, satkfunc, satn, appl_rate, ws->surfh, soil);

            infil = kinf * dh_dz;
            infil = MAX(infil, 0.0);
        }

        infil_max = appl_rate + ((ws0->surf > 0.0) ? ws0->surf / dt : 0.0);

        infil = MIN(infil, infil_max);

        infil *= wet_frac;
    }

    return infil;
}

double Recharge(const soil_struct *soil, const wstate_struct *ws, const wflux_struct *wf)
{
    double          satn;
    double          satkfunc;
    double          dh_dz;
    double          psi_u;
    double          kavg;
    double          deficit;
    double          recharge;

    if (ws->gw > soil->depth)
    {
        recharge = MIN(wf->infil, 0.0);
    }
    else if (ws->gw > soil->depth - soil->dinf && ws->unsat > 0.0)
    {
        recharge = wf->infil;
    }
    else
    {
        deficit = soil->depth - ws->gw;
        satn = ws->unsat / deficit;
        satn = MIN(satn, 1.0);
        satn = MAX(satn, SATMIN);

        satkfunc = KrFunc(soil->beta, satn);

        psi_u = Psi(satn, soil->alpha, soil->beta);

        dh_dz = (0.5 * deficit + psi_u) / (0.5 * (deficit + ws->gw));

        kavg = AvgKv(ws->gw, satkfunc, soil);

        recharge = kavg * dh_dz;

        recharge = (recharge > 0.0 && ws->unsat <= 0.0) ? 0.0 : recharge;
        recharge = (recharge < 0.0 && ws->gw <= 0.0) ? 0.0 : recharge;
    }

    return recharge;
}

double AvgKv(double gw, double satkfunc, const soil_struct *soil)
{
    double          k1, k2, k3;
    double          d1, d2, d3;
    double          deficit;

    gw = MAX(gw, 0.0);
    deficit = MAX(soil->depth - gw, 0.0);

    if (deficit > soil->dmac)
    {
        k1 = satkfunc * soil->ksatv;
        d1 = soil->dmac;

        k2 = satkfunc * soil->ksatv;
        d2 = deficit - soil->dmac;

        k3 = soil->ksatv;
        d3 = gw;
    }
    else
    {
        k1 = satkfunc * soil->ksatv;
        d1 = deficit;

        k2 = (soil->areafh > 0.0) ? soil->kmacv * soil->areafh + soil->ksatv * (1.0 - soil->areafh) : soil->ksatv;
        d2 = soil->dmac - deficit;

        k3 = soil->ksatv;
        d3 = gw - (soil->dmac - deficit);
    }

#if defined(_ARITH_)
    // Arithmetic mean formulation
    return (k1 * d1 + k2 * d2 + k3 * d3) / (d1 + d2 + d3);
#else
    return (d1 + d2 + d3) / (d1 / k1 + d2 / k2 + d3 / k3);
#endif
}

// For infiltration, macropores act as cracks, and are hydraulically effective in rapidly conducting water flow (Chen
// and Wagenet, 1992, Journal of Hydrology, 130, 105-126). For macropore hydraulic conductivities, use van Genuchten
// parameters for fractures (Gerke and van Genuchten, 1993, Water Resources Research, 29, 1225-1238).
double EffKinf(double dh_dz, double ksatfunc, double satn, double appl_rate, double surfh, const soil_struct *soil)
{
    double          keff = 0.0;
    double          kmax;
#if TEMP_DISABLED
    const double    ALPHA_CRACK = 10.0;
#endif
    const double    BETA_CRACK = 2.0;

    if (soil->areafh == 0.0)
    {
        // Matrix
        keff = soil->kinfv * ksatfunc;
    }
    else if (surfh > DEPRSTG)
    {
        // When surface wet fraction is larger than 1 (surface is totally ponded), i.e., surfh > DEPRSTG, flow situation
        // is macropore control, regardless of the application rate
        keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc + soil->kmacv * soil->areafh;
    }
    else
    {
        if (appl_rate <= dh_dz * soil->kinfv * ksatfunc)
        {
            // Matrix control
            keff = soil->kinfv * ksatfunc;
        }
        else
        {
            kmax = dh_dz * (soil->kmacv * soil->areafh + soil->kinfv * (1.0 - soil->areafh) * ksatfunc);
            if (appl_rate < kmax)
            {
                // Application control
                keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc +
                    soil->kmacv * soil->areafh * KrFunc(BETA_CRACK, satn);
            }
            else
            {
                // Macropore control
                keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc + soil->kmacv * soil->areafh;
            }
        }
    }

    return keff;
}

double Psi(double satn, double alpha, double beta)
{
    satn = MAX(satn, SATMIN);

    // van Genuchten 1980 SSSAJ
    return -pow(pow(1.0 / satn, beta / (beta - 1.0)) - 1.0, 1.0 / beta) / alpha;
}

#if defined(_DGW_)
// Hydrology for deep zone
double GeolInfil(const topo_struct *topo, const soil_struct *soil, const soil_struct *geol, const wstate_struct *ws)
{
    double          deficit;
    double          satn;
    double          psi_u;
    double          h_u;
    double          dh_dz;
    double          kavg;
    double          infil;

    if (ws->gw_geol >= geol->depth)
    {
        infil = -soil->ksatv;
    }
    else
    {
        if (ws->unsat_geol + ws->gw_geol > geol->depth || ws->gw <= 0.0)
        {
            infil = 0.0;
        }
        else
        {
            double          ksoil, kgeol;

            deficit = geol->depth - ws->gw_geol;

            satn = ws->unsat_geol / deficit;
            satn = MIN(satn, 1.0);
            satn = MAX(satn, SATMIN);

            psi_u = Psi(satn, geol->alpha, geol->beta);
            psi_u = MAX(psi_u, PSIMIN);

            h_u = psi_u + topo->zmin - 0.5 * deficit;

            dh_dz = (topo->zmin + ws->gw - h_u) / (0.5 * (ws->gw + deficit));

            ksoil = (soil->dmac >= soil->depth) ?
                soil->kmacv * soil->areafh + soil->ksatv * (1.0 - soil->areafh) : soil->ksatv;
            kgeol = (geol->dmac > 0.0) ?
                geol->ksatv * (1.0 - geol->areafh) * KrFunc(geol->beta, satn) + geol->kmacv * geol->areafh *
                KrFunc(geol->beta, satn) : geol->ksatv * KrFunc(geol->beta, satn);

            kavg = (ws->gw + deficit) / (ws->gw / ksoil + deficit / kgeol);
            infil = kavg * dh_dz;
        }
    }

    return infil;
}

double GeolRecharge(const soil_struct *geol, const wstate_struct *ws, const wflux_struct *wf)
{
    double          deficit;
    double          satn;
    double          psi_u;
    double          dh_dz;
    double          kavg;
    double          recharge;

    if (ws->gw_geol >= geol->depth)
    {
        recharge = wf->infil_geol;
    }
    else
    {
        deficit = geol->depth - ws->gw_geol;

        satn = ws->unsat_geol / deficit;
        satn = MIN(satn, 1.0);
        satn = MAX(satn, SATMIN);

        psi_u = Psi(satn, geol->alpha, geol->beta);
        psi_u = MAX(psi_u, PSIMIN);

        dh_dz = (0.5 * deficit + psi_u) / (0.5 * (deficit + ws->gw_geol));

        kavg = AvgKv(ws->gw_geol, KrFunc(geol->beta, satn), geol);

        recharge = kavg * dh_dz;

        recharge = (recharge > 0.0 && ws->unsat_geol <= 0.0) ? 0.0 : recharge;
        recharge = (recharge < 0.0 && ws->gw_geol <= 0.0) ? 0.0 : recharge;
    }

    return recharge;
}
#endif
