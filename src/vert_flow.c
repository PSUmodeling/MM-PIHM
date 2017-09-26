#include "pihm.h"

void VerticalFlow(elem_struct *elem, double dt)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          satn;
        double          satkfunc;
        double          dh_by_dz;
        double          psi_u;
        double          h_u;
        double          kinf;
        double          kavg;
        double          deficit;
        double          applrate;
        double          wetfrac;

        applrate = elem[i].wf.pcpdrp + elem[i].ws.surf / dt;

        wetfrac = (elem[i].ws.surfh > DEPRSTG) ? 1.0 :
            ((elem[i].ws.surfh < 0.0) ? 0.0 : elem[i].ws.surfh / DEPRSTG);

        if (elem[i].ws.gw > elem[i].soil.depth - elem[i].soil.dinf)
        {
            /* Assumption: Dinf < Dmac */
            dh_by_dz = (elem[i].ws.surfh + elem[i].topo.zmax -
                (elem[i].ws.gw + elem[i].topo.zmin)) / elem[i].soil.dinf;
            dh_by_dz = (elem[i].ws.surfh < 0.0 && dh_by_dz > 0.0) ?
                0.0 : dh_by_dz;
            dh_by_dz = (dh_by_dz < 1.0 && dh_by_dz > 0.0) ? 1.0 : dh_by_dz;

            satn = 1.0;
            satkfunc = KrFunc(elem[i].soil.alpha, elem[i].soil.beta, satn);

            if (elem[i].soil.areafh == 0.0)
            {
                elem[i].ps.macpore_status = MTX_CTRL;
            }
            else
            {
                elem[i].ps.macpore_status = MacroporeStatus(dh_by_dz, satkfunc,
                    applrate, elem[i].soil.kmacv, elem[i].soil.kinfv,
                    elem[i].soil.areafh);
            }

            if (dh_by_dz < 0.0)
            {
                kinf = elem[i].soil.kmacv * elem[i].soil.areafh +
                    elem[i].soil.kinfv * (1.0 - elem[i].soil.areafh);
            }
            else
            {
                kinf = EffKinf(satkfunc, satn, elem[i].ps.macpore_status,
                    elem[i].soil.kmacv, elem[i].soil.kinfv,
                    elem[i].soil.areafh);
            }

            elem[i].wf.infil = kinf * dh_by_dz;

            /* Note: infiltration can be negative in this case */
            elem[i].wf.infil = (elem[i].wf.infil < applrate) ?
                elem[i].wf.infil : applrate;

            elem[i].wf.infil *= wetfrac;

#ifdef _NOAH_
            elem[i].wf.infil *= elem[i].ps.fcr;
#endif

            elem[i].wf.rechg = elem[i].wf.infil;
        }
        else
        {
            deficit = elem[i].soil.depth - elem[i].ws.gw;
#ifdef _NOAH_
            satn = (elem[i].ws.sh2o[0] - elem[i].soil.smcmin) /
                (elem[i].soil.smcmax - elem[i].soil.smcmin);
#else
            satn = elem[i].ws.unsat / deficit;
#endif
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            psi_u = Psi(satn, elem[i].soil.alpha, elem[i].soil.beta);
            /* Note: for psi calculation using van Genuchten relation, cutting
             * the psi-sat tail at small saturation can be performed for
             * computational advantage. If you do not want to perform this,
             * comment the statement that follows */
            psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

            h_u = psi_u + elem[i].topo.zmax - 0.5 * elem[i].soil.dinf;
            dh_by_dz = (0.5 * elem[i].ws.surfh + elem[i].topo.zmax - h_u) /
                (0.5 * (elem[i].ws.surfh + elem[i].soil.dinf));
            dh_by_dz = (elem[i].ws.surfh < 0.0 && dh_by_dz > 0.0) ?
                0.0 : dh_by_dz;

            satkfunc = KrFunc(elem[i].soil.alpha, elem[i].soil.beta, satn);

            if (elem[i].soil.areafh == 0.0)
            {
                elem[i].ps.macpore_status = MTX_CTRL;
            }
            else
            {
                elem[i].ps.macpore_status = MacroporeStatus(dh_by_dz, satkfunc,
                    applrate, elem[i].soil.kmacv, elem[i].soil.kinfv,
                    elem[i].soil.areafh);
            }

            kinf = EffKinf(satkfunc, satn, elem[i].ps.macpore_status,
                elem[i].soil.kmacv, elem[i].soil.kinfv, elem[i].soil.areafh);

            elem[i].wf.infil = kinf * dh_by_dz;

            elem[i].wf.infil = (elem[i].wf.infil < applrate) ?
                elem[i].wf.infil : applrate;
            elem[i].wf.infil = (elem[i].wf.infil > 0.0) ?
                elem[i].wf.infil : 0.0;

            elem[i].wf.infil *= wetfrac;

#ifdef _NOAH_
            elem[i].wf.infil *= elem[i].ps.fcr;
#endif

            /* Arithmetic mean formulation */
            satn = elem[i].ws.unsat / deficit;
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            satkfunc = KrFunc(elem[i].soil.alpha, elem[i].soil.beta, satn);

            psi_u = Psi(satn, elem[i].soil.alpha, elem[i].soil.beta);

            dh_by_dz =
                (0.5 * deficit + psi_u) / (0.5 * (deficit + elem[i].ws.gw));

            kavg = AvgKV(elem[i].soil.dmac, deficit, elem[i].ws.gw,
                elem[i].ps.macpore_status, satkfunc, elem[i].soil.kmacv,
                elem[i].soil.ksatv, elem[i].soil.areafh);

            elem[i].wf.rechg = (deficit <= 0.0) ? 0.0 : kavg * dh_by_dz;

            elem[i].wf.rechg =
                (elem[i].wf.rechg > 0.0 && elem[i].ws.unsat <= 0.0) ?
                0.0 : elem[i].wf.rechg;
            elem[i].wf.rechg =
                (elem[i].wf.rechg < 0.0 && elem[i].ws.gw <= 0.0) ?
                0.0 : elem[i].wf.rechg;
        }
    }
}

double AvgKV(double dmac, double deficit, double gw, double macp_status,
    double satkfunc, double kmacv, double ksatv, double areafh)
{
    double          k1, k2, k3;
    double          d1, d2, d3;

    if (deficit > dmac)
    {
        k1 = EffKV(satkfunc, macp_status, kmacv, ksatv, areafh);
        d1 = dmac;

        k2 = satkfunc * ksatv;
        d2 = deficit - dmac;

        k3 = ksatv;
        d3 = gw;
    }
    else
    {
        k1 = EffKV(satkfunc, macp_status, kmacv, ksatv, areafh);
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

double EffKV(double ksatfunc, int status, double mackv, double kv, double areaf)
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
