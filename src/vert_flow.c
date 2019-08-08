#include "pihm.h"

void VerticalFlow(elem_struct *elem, double dt)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Calculate infiltration rate */
        elem[i].wf.infil = Infil(&elem[i].ws, &elem[i].ws0, &elem[i].wf,
            &elem[i].topo, &elem[i].soil, dt);

#if defined(_NOAH_)
        /* Urban surface is assumed to be impermeable */
        elem[i].wf.infil *= 1.0 - elem[i].lc.urban_fraction;

        /* Constrain infiltration by frozen top soil */
        elem[i].wf.infil *= elem[i].ps.fcr;
#endif

        /* Calculate recharge rate */
        elem[i].wf.rechg = Recharge(&elem[i].ws, &elem[i].wf, &elem[i].soil);

#if defined(_FBR_)
        elem[i].wf.fbr_infil = FbrInfil(&elem[i].ws, &elem[i].soil,
            &elem[i].geol, &elem[i].topo);
        elem[i].wf.fbr_rechg = FbrRecharge(&elem[i].ws, &elem[i].wf,
            &elem[i].geol);
#endif
    }
}

double Infil(const wstate_struct *ws, const wstate_struct *ws0,
    const wflux_struct *wf, const topo_struct *topo, const soil_struct *soil,
    double dt)
{
    double          applrate;
    double          wetfrac;
    double          dh_by_dz;
    double          satn;
    double          satkfunc;
    double          infil;
    double          infil_max;
    double          kinf;
    double          deficit;
    double          psi_u;
    double          h_u;
    int             j;

    if (ws->unsat + ws->gw > soil->depth)
    {
        infil = 0.0;
    }
    else
    {
        /* Water appliation rate is the sum of net gain of overland flow and net
         * precipitation */
        applrate = 0.0;
        for (j = 0; j < NUM_EDGE; j++)
        {
            applrate += -wf->ovlflow[j] / topo->area;
        }
        applrate = (applrate > 0.0) ? applrate : 0.0;
        applrate += wf->pcpdrp;

        if (DEPRSTG > 0.0)
        {
            wetfrac = ws->surfh / DEPRSTG;
            wetfrac = (wetfrac > 0.0) ? wetfrac : 0.0;
            wetfrac = (wetfrac < 1.0) ? wetfrac : 1.0;
        }
        else
        {
            wetfrac = (ws->surfh > 0.0) ? 1.0 : 0.0;
        }

        if (ws->gw > soil->depth - soil->dinf)
        {
            /* Assumption: Dinf < Dmac */
            dh_by_dz = (ws->surfh + topo->zmax - (ws->gw + topo->zmin)) /
                (0.5 * (ws->surfh + soil->dinf));
            dh_by_dz = (ws->surfh <= 0.0 && dh_by_dz > 0.0) ? 0.0 : dh_by_dz;

            satn = 1.0;
            satkfunc = KrFunc(soil->beta, satn);

            kinf = EffKinf(soil, dh_by_dz, satkfunc, satn, applrate, ws->surfh);

            infil = kinf * dh_by_dz;
        }
        else
        {
            deficit = soil->depth - ws->gw;
#if defined(_NOAH_)
            satn = (ws->sh2o[0] - soil->smcmin) / (soil->smcmax - soil->smcmin);
#else
            satn = ws->unsat / deficit;
#endif
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            psi_u = Psi(satn, soil->alpha, soil->beta);
            /* Note: for psi calculation using van Genuchten relation, cutting
             * the psi-sat tail at small saturation can be performed for
             * computational advantage. If you do not want to perform this,
             * comment the statement that follows */
            psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

            h_u = psi_u + topo->zmax - 0.5 * soil->dinf;
            dh_by_dz = (ws->surfh + topo->zmax - h_u) /
                (0.5 * (ws->surfh + soil->dinf));
            dh_by_dz = (ws->surfh <= 0.0 && dh_by_dz > 0.0) ?  0.0 : dh_by_dz;

            satkfunc = KrFunc(soil->beta, satn);

            kinf = EffKinf(soil, dh_by_dz, satkfunc, satn, applrate, ws->surfh);

            infil = kinf * dh_by_dz;
            infil = (infil > 0.0) ? infil : 0.0;
        }

        infil_max = applrate + ((ws0->surf > 0.0) ? ws0->surf / dt : 0.0);

        infil = (infil > infil_max) ? infil_max : infil;

        infil *= wetfrac;
    }

    return infil;
}

double Recharge(const wstate_struct *ws, const wflux_struct *wf,
    const soil_struct *soil)
{
    double          satn;
    double          satkfunc;
    double          dh_by_dz;
    double          psi_u;
    double          kavg;
    double          deficit;
    double          rechg;

    if (ws->gw > soil->depth - soil->dinf)
    {
        rechg = wf->infil;
    }
    else
    {
        deficit = soil->depth - ws->gw;
        satn = ws->unsat / deficit;
        satn = (satn > 1.0) ? 1.0 : satn;
        satn = (satn < SATMIN) ? SATMIN : satn;

        satkfunc = KrFunc(soil->beta, satn);

        psi_u = Psi(satn, soil->alpha, soil->beta);

        dh_by_dz =
            (0.5 * deficit + psi_u) / (0.5 * (deficit + ws->gw));

        kavg = AvgKv(soil, deficit, ws->gw, satkfunc);

        rechg = kavg * dh_by_dz;

        rechg = (rechg > 0.0 && ws->unsat <= 0.0) ?  0.0 : rechg;
        rechg = (rechg < 0.0 && ws->gw <= 0.0) ?  0.0 : rechg;
    }

    return rechg;
}

double AvgKv(const soil_struct *soil, double deficit, double gw,
    double satkfunc)
{
    double          k1, k2, k3;
    double          d1, d2, d3;

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

        k2 = (soil->areafh > 0.0) ?
            soil->kmacv * soil->areafh + soil->ksatv * (1.0 - soil->areafh) :
            soil->ksatv;
        d2 = soil->dmac - deficit;

        k3 = soil->ksatv;
        d3 = gw - (soil->dmac - deficit);
    }

#if defined(_ARITH_)
    /* Arithmetic mean formulation */
    return (k1 * d1 + k2 * d2 + k3 * d3) / (d1 + d2 + d3);
#else
    return (d1 + d2 + d3) / (d1 / k1 + d2 / k2 + d3 / k3);
#endif
}

double EffKinf(const soil_struct *soil, double dh_by_dz, double ksatfunc,
    double elemsatn, double applrate, double surfh)
{
    /*
     * For infiltration, macropores act as cracks, and are hydraulically
     * effective in rapidly conducting water flow (Chen and Wagenet, 1992,
     * Journal of Hydrology, 130, 105-126).
     * For macropore hydraulic conductivities, use van Genuchten parameters for
     * fractures (Gerke and van Genuchten, 1993, Water Resources Research, 29,
     * 1225-1238).
     */
    double          keff = 0.0;
    double          kmax;
#if TEMP_DISABLED
    const double    ALPHA_CRACK = 10.0;
#endif
    const double    BETA_CRACK = 2.0;

    if (soil->areafh == 0.0)
    {
        /* Matrix */
        keff = soil->kinfv * ksatfunc;
    }
    else if (surfh > DEPRSTG)
    {
        /* When surface wet fraction is larger than 1 (surface is totally
         * ponded), i.e., surfh > DEPRSTG, flow situation is macropore control,
         * regardless of the application rate */
        keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc +
            soil->kmacv * soil->areafh;
    }
    else
    {
        if (applrate <= dh_by_dz * soil->kinfv * ksatfunc)
        {
            /* Matrix control */
            keff = soil->kinfv * ksatfunc;
        }
        else
        {
            kmax = dh_by_dz * (soil->kmacv * soil->areafh +
                soil->kinfv * (1.0 - soil->areafh) * ksatfunc);
            if (applrate < kmax)
            {
                /* Application control */
                keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc +
                    soil->kmacv * soil->areafh *
                    KrFunc(BETA_CRACK, elemsatn);
            }
            else
            {
                /* Macropore control */
                keff = soil->kinfv * (1.0 - soil->areafh) * ksatfunc +
                    soil->kmacv * soil->areafh;
            }
        }
    }

    return keff;
}

double Psi(double satn, double alpha, double beta)
{
    satn = (satn < SATMIN) ? SATMIN : satn;

    /* van Genuchten 1980 SSSAJ */
    return -pow(pow(1.0 / satn, beta / (beta - 1.0)) - 1.0, 1.0 / beta) / alpha;
}

#if defined(_FBR_)
/*
 * Hydrology for fractured bedrock
 */
double FbrInfil(const wstate_struct *ws, const soil_struct *soil,
    const geol_struct *geol, const topo_struct *topo)
{
    double          deficit;
    double          satn;
    double          psi_u;
    double          h_u;
    double          satkfunc;
    double          dh_by_dz;
    double          kavg;
    double          infil;

    if (ws->fbr_gw >= geol->depth)
    {
# if defined(_TGM_)
        /* In two-grid model, oversaturated deep groundwater should enter stream
         * instead of shallow groundwater */
        infil = 0.0;
# else
        infil = -soil->ksatv;
# endif
    }
    else
    {
        if (ws->fbr_unsat + ws->fbr_gw > geol->depth || ws->gw <= 0.0)
        {
            infil = 0.0;
        }
        else
        {
            deficit = geol->depth - ws->fbr_gw;

            satn = ws->fbr_unsat / deficit;
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            psi_u = Psi(satn, geol->alpha, geol->beta);
            psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

            h_u = psi_u + topo->zmin - 0.5 * deficit;

            satkfunc = KrFunc(geol->beta, satn);

            dh_by_dz = (topo->zmin + ws->gw - h_u) / (0.5 * (ws->gw + deficit));

            kavg = (ws->gw + deficit) /
                (ws->gw / soil->ksatv + deficit / (geol->ksatv * satkfunc));
            infil = kavg * dh_by_dz;
        }
    }

    return infil;
}

double FbrRecharge(const wstate_struct *ws, const wflux_struct *wf,
    const geol_struct *geol)
{
    double          deficit;
    double          satn;
    double          psi_u;
    double          satkfunc;
    double          dh_by_dz;
    double          kavg;
    double          rechg;

    if (ws->fbr_gw >= geol->depth)
    {
# if defined(_TGM_)
        rechg = 0.0;
# else
        rechg = wf->fbr_infil;
# endif
    }
    else
    {
        deficit = geol->depth - ws->fbr_gw;

        satn = ws->fbr_unsat / deficit;
        satn = (satn > 1.0) ? 1.0 : satn;
        satn = (satn < SATMIN) ? SATMIN : satn;

        psi_u = Psi(satn, geol->alpha, geol->beta);
        psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

        satkfunc = KrFunc(geol->beta, satn);

        dh_by_dz = (0.5 * deficit + psi_u) / (0.5 * (deficit + ws->fbr_gw));

        kavg = (ws->fbr_unsat * geol->ksatv * satkfunc +
             ws->fbr_gw * geol->ksatv) /
            (ws->fbr_unsat + ws->fbr_gw);

        rechg = kavg * dh_by_dz;

        rechg = (rechg > 0.0 && ws->fbr_unsat <= 0.0) ? 0.0 : rechg;
        rechg = (rechg < 0.0 && ws->fbr_gw <= 0.0) ? 0.0 : rechg;
    }

    return rechg;
}
#endif
