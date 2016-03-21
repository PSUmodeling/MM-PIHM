#include "pihm.h"

void VerticalFlow (pihm_struct pihm)
{
    int             i;
    double          satn;
    double          satkfunc;
    double          dh_by_dz;
    double          psi_u;
    double          h_u;
    double          kinf;
    double          kavg;
    double          dt;
    double          deficit;
    double          applrate;
    double          wetfrac;
    elem_struct    *elem;

    dt = (double)pihm->ctrl.stepsize;

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        applrate = elem->wf.netprcp + elem->ws.surf / dt;

        wetfrac =
            (elem->ws.surf > DEPRSTG) ? 1.0 : ((elem->ws.surf <
                0.0) ? 0.0 : pow (elem->ws.surf / DEPRSTG, 2.0));

        if (elem->ws.gw > elem->soil.depth - elem->soil.dinf)
        {
            /* Assumption: Dinf < Dmac */
            dh_by_dz =
                (elem->ws.surf + elem->topo.zmax - (elem->ws.gw +
                    elem->topo.zmin)) / elem->soil.dinf;
            dh_by_dz = (elem->ws.surf < 0.0 &&
                dh_by_dz > 0.0) ? 0.0 : dh_by_dz;
            dh_by_dz = (dh_by_dz < 1.0 && dh_by_dz > 0.0) ? 1.0 : dh_by_dz;

            satn = 1.0;
            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            if (elem->soil.macropore)
            {
                elem->ps.macpore_status =
                    MacroporeStatus (dh_by_dz, satkfunc, satn, applrate,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }
            else
            {
                elem->ps.macpore_status = MTX_CTRL;
            }

            if (dh_by_dz < 0.0)
            {
                kinf = (elem->soil.macropore) ?
                    elem->soil.kmacv * elem->soil.areafh +
                    elem->soil.kinfv * (1.0 - elem->soil.areafh) :
                    elem->soil.kinfv;
            }
            else
            {
                kinf = EffKinf (satkfunc, satn, elem->ps.macpore_status,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }

            elem->wf.infil = kinf * dh_by_dz;

            /* Note: infiltration can be negative in this case */
            elem->wf.infil =
                (elem->wf.infil < applrate) ? elem->wf.infil : applrate;

            elem->wf.infil *= wetfrac;

#ifdef _NOAH_
            elem->wf.infil *= elem->ps.fcr;
#endif


            elem->wf.rechg = elem->wf.infil;
        }
        else
        {
            deficit = elem->soil.depth - elem->ws.gw;
#ifdef _NOAH_
            satn = (elem->ws.sh2o[0] - elem->soil.smcmin) /
                (elem->soil.smcmax - elem->soil.smcmin);
#else
            satn = elem->ws.unsat / deficit;
#endif
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;
            /* Note: for psi calculation using van genuchten relation, cutting
             * the psi-sat tail at small saturation can be performed for
             * computational advantage. if you dont' want to perform this,
             * comment the statement that follows */
            psi_u = Psi (satn, elem->soil.alpha, elem->soil.beta);
            psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

            h_u = psi_u + elem->topo.zmax - 0.5 * elem->soil.dinf;
            dh_by_dz =
                (0.5 * elem->ws.surf + elem->topo.zmax -
                h_u) / (0.5 * (elem->ws.surf + elem->soil.dinf));
            dh_by_dz = (elem->ws.surf < 0.0 &&
                dh_by_dz > 0.0) ? 0.0 : dh_by_dz;

            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            if (elem->soil.macropore)
            {
                elem->ps.macpore_status =
                    MacroporeStatus (dh_by_dz, satkfunc, satn, applrate,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }
            else
            {
                elem->ps.macpore_status = MTX_CTRL;
            }

            kinf = EffKinf (satkfunc, satn, elem->ps.macpore_status,
                elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);

            elem->wf.infil = kinf * dh_by_dz;

            elem->wf.infil =
                (elem->wf.infil < applrate) ? elem->wf.infil : applrate;
            elem->wf.infil = (elem->wf.infil > 0.0) ? elem->wf.infil : 0.0;

            elem->wf.infil *= wetfrac;

#ifdef _NOAH_
            elem->wf.infil *= elem->ps.fcr;
#endif

            /* Harmonic mean formulation.
             * Note that if unsaturated zone has low saturation, satkfunc
             * becomes very small. use arithmetic mean instead */
            //elem->rechg = (satn==0.0)?0:(deficit<=0)?0:(md->ele[i].ksatv*satkfunc*(md->ele[i].alpha*deficit-2*pow(-1+pow(satn,md->ele[i].beta/(-md->ele[i].beta+1)),1/md->ele[i].beta))/(md->ele[i].alpha*((deficit+md->dummyy[i+2*md->numele]*satkfunc))));
            /* Arithmetic mean formulation */
            satn = elem->ws.unsat / deficit;
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            //if (elem->ps.macpore_status > MTX_CTRL && elem->ws.gw > elem->soil.depth - elem->soil.dmac)
            //{
            //    keff = EffKV (satkfunc, satn, elem->ps.macpore_status, elem->soil.kmacv, elem->soil.ksatv, elem->soil.areafh);
            //}
            //else
            //{
            //    keff = elem->soil.ksatv * satkfunc;
            //}

            psi_u = Psi (satn, elem->soil.alpha, elem->soil.beta);

            dh_by_dz =
                (0.5 * deficit + psi_u) / (0.5 * (deficit + elem->ws.gw));

            //kavg = 1.0 / (deficit / effk[KMTX] + elem->ws.gw / elem->soil.ksatv);
            //kavg = (deficit * keff + elem->ws.gw * elem->soil.ksatv) / (deficit + elem->ws.gw);
            kavg =
                AvgKV (elem->soil.dmac, deficit, elem->ws.gw,
                elem->ps.macpore_status, satn, satkfunc, elem->soil.kmacv,
                elem->soil.ksatv, elem->soil.areafh);

            elem->wf.rechg = (deficit <= 0.0) ? 0.0 : kavg * dh_by_dz;

            elem->wf.rechg = (elem->wf.rechg > 0.0 &&
                elem->ws.unsat <= 0.0) ? 0.0 : elem->wf.rechg;
            elem->wf.rechg = (elem->wf.rechg < 0.0 &&
                elem->ws.gw <= 0.0) ? 0.0 : elem->wf.rechg;
        }
    }
}

double AvgKV (double dmac, double deficit, double gw, double macp_status,
    double satn, double satkfunc, double kmacv, double ksatv, double areafh)
{
    double          k1, k2, k3;
    double          d1, d2, d3;

    if (deficit > dmac)
    {
        k1 = EffKV (satkfunc, satn, macp_status, kmacv, ksatv, areafh);
        d1 = dmac;

        k2 = satkfunc * ksatv;
        d2 = deficit - dmac;

        k3 = ksatv;
        d3 = gw;
    }
    else
    {
        k1 = EffKV (satkfunc, satn, macp_status, kmacv, ksatv, areafh);
        d1 = deficit;

        k2 = kmacv * areafh + ksatv * (1.0 - areafh);
        d2 = dmac - deficit;

        k3 = ksatv;
        d3 = gw - (dmac - deficit);
    }

#ifdef _ARITH_
    return ((k1 * d1 + k2 * d2 + k3 * d3) / (d1 + d2 + d3));
#else
    return (d1 + d2 + d3) / (d1 / k1 + d2 / k2 + d3 / k3);
#endif
}

double EffKinf (double ksatfunc, double elemsatn, int status, double mackv,
    double kinf, double areaf)
{
    double          keff = 0.0;

    switch (status)
    {
        case MTX_CTRL:
            keff = kinf * ksatfunc;
            break;
        case APP_CTRL:
            keff =
                kinf * (1.0 - areaf) * ksatfunc +
                mackv * areaf * KrFunc (10.0, 2.0, elemsatn);
            break;
        case MAC_CTRL:
            keff = kinf * (1.0 - areaf) * ksatfunc + mackv * areaf;
            break;
        default:
            printf ("Error: Macropore status not recognized!\n");
            PihmExit (1);
    }

    return (keff);
}

double EffKV (double ksatfunc, double elemsatn, int status, double mackv,
    double kv, double areaf)
{
    double          keff = 0.0;

    switch (status)
    {
        case MTX_CTRL:
            keff = kv * ksatfunc;
            break;
        case APP_CTRL:
            keff = kv * (1.0 - areaf) * ksatfunc + mackv * areaf * ksatfunc;    //KrFunc (10.0, 2.0, elemsatn);
            break;
        case MAC_CTRL:
            keff = kv * (1.0 - areaf) * ksatfunc + mackv * areaf;
            break;
        default:
            printf ("Error: Macropore status not recognized!\n");
            PihmExit (1);
    }

    return (keff);
}

double KrFunc (double alpha, double beta, double satn)
{
    return (pow (satn, 0.5) * pow (-1.0 + pow (1.0 - pow (satn,
                    beta / (beta - 1.0)), (beta - 1.0) / beta), 2));
}

int MacroporeStatus (double dh_by_dz, double ksatfunc, double elemsatn,
    double applrate, double mackv, double kv, double areaf)
{
    dh_by_dz = (dh_by_dz < 1.0) ? 1.0 : dh_by_dz;

    if (applrate <= dh_by_dz * kv * ksatfunc)
    {
        return (MTX_CTRL);
    }
    else
    {
        if (applrate <
            dh_by_dz * (mackv * areaf + kv * (1.0 - areaf) * ksatfunc))
        {
            return (APP_CTRL);
        }
        else
        {
            return (MAC_CTRL);
        }
    }
}

double Psi (double satn, double alpha, double beta)
{
    /* van Genuchten 1980 SSSAJ */
    return (0.0 - pow (pow (1.0 / satn, beta / (beta - 1.0)) - 1.0,
            1.0 / beta) / alpha);
}
