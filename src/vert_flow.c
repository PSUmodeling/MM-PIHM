#include "pihm.h"

void VerticalFlow(elem_struct *elem, double dt)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].wf.infil = Infil(&elem[i].ws, &elem[i].wf, &elem[i].topo,
            &elem[i].soil, dt, &elem[i].ps);

        elem[i].wf.rechg = Recharge(&elem[i].ws, &elem[i].wf, &elem[i].ps,
            &elem[i].soil);
    }
}

double Infil(const wstate_struct *ws, const wflux_struct *wf,
    const topo_struct *topo, const soil_struct *soil, double dt,
    pstate_struct *ps)
{
    double          applrate;
    double          wetfrac;
    double          dh_by_dz;
    double          satn;
    double          satkfunc;
    double          infil;
    double          kinf;
    double          deficit;
    double          psi_u;
    double          h_u;

    applrate = wf->pcpdrp + ws->surf / dt;

    wetfrac = ws->surfh / DEPRSTG;
    wetfrac = (wetfrac > 0.0) ? wetfrac : 0.0;
    wetfrac = (wetfrac < 1.0) ? wetfrac : 1.0;

    if (ws->gw > soil->depth - soil->dinf)
    {
        /* Assumption: Dinf < Dmac */
        dh_by_dz =
            (ws->surfh + topo->zmax - (ws->gw + topo->zmin)) / soil->dinf;
        dh_by_dz = (ws->surfh < 0.0 && dh_by_dz > 0.0) ? 0.0 : dh_by_dz;
        dh_by_dz = (dh_by_dz < 1.0 && dh_by_dz > 0.0) ? 1.0 : dh_by_dz;

        satn = 1.0;
        satkfunc = KrFunc(soil->alpha, soil->beta, satn);

        if (soil->areafh == 0.0)
        {
            ps->macpore_status = MTX_CTRL;
        }
        else
        {
            ps->macpore_status = MacroporeStatus(dh_by_dz, satkfunc,
                applrate, soil->kmacv, soil->kinfv, soil->areafh);
        }

        if (dh_by_dz < 0.0)
        {
            kinf = soil->kmacv * soil->areafh +
                soil->kinfv * (1.0 - soil->areafh);
        }
        else
        {
            kinf = EffKinf(satkfunc, satn, ps->macpore_status,
                soil->kmacv, soil->kinfv, soil->areafh);
        }

        infil = kinf * dh_by_dz;
    }
    else
    {
        deficit = soil->depth - ws->gw;
#ifdef _NOAH_
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
        dh_by_dz = (0.5 * ws->surfh + topo->zmax - h_u) /
            (0.5 * (ws->surfh + soil->dinf));
        dh_by_dz = (ws->surfh < 0.0 && dh_by_dz > 0.0) ?  0.0 : dh_by_dz;

        satkfunc = KrFunc(soil->alpha, soil->beta, satn);

        if (soil->areafh == 0.0)
        {
            ps->macpore_status = MTX_CTRL;
        }
        else
        {
            ps->macpore_status = MacroporeStatus(dh_by_dz, satkfunc,
                applrate, soil->kmacv, soil->kinfv, soil->areafh);
        }

        kinf = EffKinf(satkfunc, satn, ps->macpore_status,
            soil->kmacv, soil->kinfv, soil->areafh);

        infil = kinf * dh_by_dz;

        infil = (infil > 0.0) ? infil : 0.0;
    }

    infil = (infil < applrate) ? infil : applrate;

    infil *= wetfrac;

#ifdef _NOAH_
    infil *= ps->fcr;
#endif

    return infil;
}

double Recharge(const wstate_struct *ws, const wflux_struct *wf,
    const pstate_struct *ps, const soil_struct *soil)
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
        /* Arithmetic mean formulation */
        deficit = soil->depth - ws->gw;
        satn = ws->unsat / deficit;
        satn = (satn > 1.0) ? 1.0 : satn;
        satn = (satn < SATMIN) ? SATMIN : satn;

        satkfunc = KrFunc(soil->alpha, soil->beta, satn);

        psi_u = Psi(satn, soil->alpha, soil->beta);

        dh_by_dz =
            (0.5 * deficit + psi_u) / (0.5 * (deficit + ws->gw));

        kavg = AvgKv(soil->dmac, deficit, ws->gw,
            ps->macpore_status, satkfunc, soil->kmacv,
            soil->ksatv, soil->areafh);

        rechg = kavg * dh_by_dz;

        rechg = (rechg > 0.0 && ws->unsat <= 0.0) ?  0.0 : rechg;
        rechg = (rechg < 0.0 && ws->gw <= 0.0) ?  0.0 : rechg;
    }

    return rechg;
}

double AvgKv(double dmac, double deficit, double gw, double macp_status,
    double satkfunc, double kmacv, double ksatv, double areafh)
{
    double          k1, k2, k3;
    double          d1, d2, d3;

    if (deficit > dmac)
    {
        k1 = EffKv(satkfunc, macp_status, kmacv, ksatv, areafh);
        d1 = dmac;

        k2 = satkfunc * ksatv;
        d2 = deficit - dmac;

        k3 = ksatv;
        d3 = gw;
    }
    else
    {
        k1 = EffKv(satkfunc, macp_status, kmacv, ksatv, areafh);
        d1 = deficit;

        k2 = kmacv * areafh + ksatv * (1.0 - areafh);
        d2 = dmac - deficit;

        k3 = ksatv;
        d3 = gw - (dmac - deficit);
    }

#ifdef _ARITH_
    return (k1 * d1 + k2 * d2 + k3 * d3) / (d1 + d2 + d3);
#else
    return (d1 + d2 + d3) / (d1 / k1 + d2 / k2 + d3 / k3);
#endif
}

double EffKinf(double ksatfunc, double elemsatn, int status, double mackv,
    double kinf, double areaf)
{
    double          keff = 0.0;
    const double    ALPHA_CRACK = 10.0;
    const double    BETA_CRACK = 2.0;

    switch (status)
    {
        case MTX_CTRL:
            keff = kinf * ksatfunc;
            break;
        case APP_CTRL:
            keff = kinf * (1.0 - areaf) * ksatfunc +
                mackv * areaf * KrFunc(ALPHA_CRACK, BETA_CRACK, elemsatn);
            break;
        case MAC_CTRL:
            keff = kinf * (1.0 - areaf) * ksatfunc + mackv * areaf;
            break;
        default:
            PIHMprintf(VL_ERROR,
                "Error: Macropore status (%d) is not defined.\n", status);
            PIHMexit(EXIT_FAILURE);
    }

    return keff;
}

double EffKv(double ksatfunc, int status, double mackv, double kv, double areaf)
{
    double          keff = 0.0;

    switch (status)
    {
        case MTX_CTRL:
            keff = kv * ksatfunc;
            break;
        case APP_CTRL:
            keff = kv * (1.0 - areaf) * ksatfunc + mackv * areaf * ksatfunc;
            break;
        case MAC_CTRL:
            keff = kv * (1.0 - areaf) * ksatfunc + mackv * areaf;
            break;
        default:
            PIHMprintf(VL_ERROR,
                "Error: Macropore status (%d) is not defined.\n", status);
            PIHMexit(EXIT_FAILURE);
    }

    return keff;
}

double KrFunc(double alpha, double beta, double satn)
{
    return sqrt(satn) *
        (1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta)) *
        (1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta));
}

int MacroporeStatus(double dh_by_dz, double ksatfunc, double applrate,
    double mackv, double kv, double areaf)
{
    dh_by_dz = (dh_by_dz < 1.0) ? 1.0 : dh_by_dz;

    if (applrate <= dh_by_dz * kv * ksatfunc)
    {
        return MTX_CTRL;
    }
    else
    {
        if (applrate <
            dh_by_dz * (mackv * areaf + kv * (1.0 - areaf) * ksatfunc))
        {
            return APP_CTRL;
        }
        else
        {
            return MAC_CTRL;
        }
    }
}

double Psi(double satn, double alpha, double beta)
{
    /* van Genuchten 1980 SSSAJ */
    return -pow(pow(1.0 / satn, beta / (beta - 1.0)) - 1.0, 1.0 / beta) / alpha;
}
