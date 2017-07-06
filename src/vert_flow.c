#include "pihm.h"

void VerticalFlow (pihm_struct pihm)
{
    int             i;
    double          dt;

    dt = (double)pihm->ctrl.stepsize;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double      satn;
        double      satkfunc;
        double      dh_by_dz;
        double      psi_u;
        double      h_u;
        double      kinf;
        double      kavg;
        double      deficit;
        double      applrate;
        double      wetfrac;
        elem_struct *elem;

        elem = &pihm->elem[i];

        applrate = elem->wf.pcpdrp + elem->ws.surf / dt;

        wetfrac = (elem->ws.surfh > DEPRSTG) ?  1.0 :
            ((elem->ws.surfh < 0.0) ?  0.0 : elem->ws.surfh / DEPRSTG);

        if (elem->ws.gw > elem->soil.depth - elem->soil.dinf)
        {
            /* Assumption: Dinf < Dmac */
            dh_by_dz = (elem->ws.surfh + elem->topo.zmax -
                (elem->ws.gw + elem->topo.zmin)) / elem->soil.dinf;
            dh_by_dz = (elem->ws.surfh < 0.0 && dh_by_dz > 0.0) ?
                0.0 : dh_by_dz;
            dh_by_dz = (dh_by_dz < 1.0 && dh_by_dz > 0.0) ? 1.0 : dh_by_dz;

            satn = 1.0;
            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            if (elem->soil.areafh == 0.0)
            {
                elem->ps.macpore_status = MTX_CTRL;
            }
            else
            {
                elem->ps.macpore_status =
                    MacroporeStatus (dh_by_dz, satkfunc, applrate,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }

            if (dh_by_dz < 0.0)
            {
                kinf = elem->soil.kmacv * elem->soil.areafh +
                    elem->soil.kinfv * (1.0 - elem->soil.areafh);
            }
            else
            {
                kinf = EffKinf (satkfunc, satn, elem->ps.macpore_status,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }

            elem->wf.infil = kinf * dh_by_dz;

            /* Note: infiltration can be negative in this case */
            elem->wf.infil = (elem->wf.infil < applrate) ?
                elem->wf.infil : applrate;

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

            psi_u = Psi (satn, elem->soil.alpha, elem->soil.beta);
            /* Note: for psi calculation using van genuchten relation, cutting
             * the psi-sat tail at small saturation can be performed for
             * computational advantage. if you dont' want to perform this,
             * comment the statement that follows */
            psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

            h_u = psi_u + elem->topo.zmax - 0.5 * elem->soil.dinf;
            dh_by_dz =
                (0.5 * elem->ws.surfh + elem->topo.zmax -
                h_u) / (0.5 * (elem->ws.surfh + elem->soil.dinf));
            dh_by_dz = (elem->ws.surfh < 0.0 && dh_by_dz > 0.0) ?
                0.0 : dh_by_dz;

            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            if (elem->soil.areafh == 0.0)
            {
                elem->ps.macpore_status = MTX_CTRL;
            }
            else
            {
                elem->ps.macpore_status =
                    MacroporeStatus (dh_by_dz, satkfunc, applrate,
                    elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);
            }

            kinf = EffKinf (satkfunc, satn, elem->ps.macpore_status,
                elem->soil.kmacv, elem->soil.kinfv, elem->soil.areafh);

            elem->wf.infil = kinf * dh_by_dz;

            elem->wf.infil = (elem->wf.infil < applrate) ?
                elem->wf.infil : applrate;
            elem->wf.infil = (elem->wf.infil > 0.0) ? elem->wf.infil : 0.0;

            elem->wf.infil *= wetfrac;

#ifdef _NOAH_
            elem->wf.infil *= elem->ps.fcr;
#endif

            /* Arithmetic mean formulation */
            satn = elem->ws.unsat / deficit;
            satn = (satn > 1.0) ? 1.0 : satn;
            satn = (satn < SATMIN) ? SATMIN : satn;

            satkfunc = KrFunc (elem->soil.alpha, elem->soil.beta, satn);

            psi_u = Psi (satn, elem->soil.alpha, elem->soil.beta);

            dh_by_dz =
                (0.5 * deficit + psi_u) / (0.5 * (deficit + elem->ws.gw));

            kavg = AvgKV (elem->soil.dmac, deficit, elem->ws.gw,
                elem->ps.macpore_status, satkfunc, elem->soil.kmacv,
                elem->soil.ksatv, elem->soil.areafh);

            elem->wf.rechg = (deficit <= 0.0) ? 0.0 : kavg * dh_by_dz;

            elem->wf.rechg = (elem->wf.rechg > 0.0 && elem->ws.unsat <= 0.0) ?
                0.0 : elem->wf.rechg;
            elem->wf.rechg = (elem->wf.rechg < 0.0 && elem->ws.gw <= 0.0) ?
                0.0 : elem->wf.rechg;
        }
    }
}

double AvgKV (double dmac, double deficit, double gw, double macp_status,
    double satkfunc, double kmacv, double ksatv, double areafh)
{
    double          k1, k2, k3;
    double          d1, d2, d3;

    if (deficit > dmac)
    {
        k1 = EffKV (satkfunc, macp_status, kmacv, ksatv, areafh);
        d1 = dmac;

        k2 = satkfunc * ksatv;
        d2 = deficit - dmac;

        k3 = ksatv;
        d3 = gw;
    }
    else
    {
        k1 = EffKV (satkfunc, macp_status, kmacv, ksatv, areafh);
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
    const double    ALPHA_CRACK = 10.0;
    const double    BETA_CRACK = 2.0;

    switch (status)
    {
        case MTX_CTRL:
            keff = kinf * ksatfunc;
            break;
        case APP_CTRL:
            keff = kinf * (1.0 - areaf) * ksatfunc +
                mackv * areaf * KrFunc (ALPHA_CRACK, BETA_CRACK, elemsatn);
            break;
        case MAC_CTRL:
            keff = kinf * (1.0 - areaf) * ksatfunc + mackv * areaf;
            break;
        default:
            PIHMprintf (VL_ERROR,
                "Error: Macropore status (%d) is not defined.\n", status);
            PIHMexit (EXIT_FAILURE);
    }

    return (keff);
}

double EffKV (double ksatfunc, int status, double mackv, double kv,
    double areaf)
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
            PIHMprintf (VL_ERROR,
                "Error: Macropore status (%d) is not defined.\n", status);
            PIHMexit (EXIT_FAILURE);
    }

    return (keff);
}

double KrFunc (double alpha, double beta, double satn)
{
    // If satn < 0.4 kr~0
	if (satn < 0.4)
	{
		return 0.0;
	}
	else
	{
		double n = (beta - 1.0) / beta;
		return (sqrt(satn) * (-1.0 + pow(1.0 - pow(satn,
			1.0 / n), n))*(-1.0 + pow(1.0 - pow(satn,
				1.0 / n), n)));
	}
}

int MacroporeStatus (double dh_by_dz, double ksatfunc, double applrate,
    double mackv, double kv, double areaf)
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
	double n1, n2, invsat, rv=0.0;

    //TODO remove this?
	if (satn < 0.9999999) 
	{ 
		n1 = beta / (beta - 1.0);
		n2 = 1.0 / beta;
		invsat = 1.0 / satn;	
		rv = -pow(pow(invsat, n1) - 1.0, n2) / alpha; 
	}
    /* van Genuchten 1980 SSSAJ */
    return (0.0 -
        pow (pow (1.0 / satn, beta / (beta - 1.0)) - 1.0,
            1.0 / beta) / alpha);
}
