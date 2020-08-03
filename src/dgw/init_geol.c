#include "pihm.h"

void InitGeol(const geoltbl_struct *geoltbl, const calib_struct *calib,
    elem_struct elem[])
{
    int             i;
    int             geol_ind;

    for (i = 0; i < nelem; i++)
    {
        if (elem[i].attrib.geol_type > geoltbl->number)
        {
            pihm_printf(VL_ERROR,
                "Error: Geol type %d for Element %d is not in the geol file.",
                elem[i].attrib.geol_type, i + 1);
            pihm_exit(EXIT_FAILURE);
        }

        geol_ind = elem[i].attrib.geol_type - 1;

        elem[i].geol.depth = elem[i].topo.zmin - elem[i].topo.zbed;

        elem[i].geol.ksath = calib->geol_ksath * geoltbl->ksath[geol_ind];
        elem[i].geol.ksatv = calib->geol_ksatv * geoltbl->ksatv[geol_ind];

        elem[i].geol.smcmin = calib->geol_porosity * geoltbl->smcmin[geol_ind];
        elem[i].geol.smcmax = calib->geol_porosity * geoltbl->smcmax[geol_ind];
        elem[i].geol.porosity = elem[i].geol.smcmax - elem[i].geol.smcmin;
        if (elem[i].geol.porosity > 1.0 || elem[i].geol.porosity <= 0.0)
        {
            pihm_printf (VL_ERROR,
                "Error: Porosity value out of bounds for Element %d", i + 1);
            pihm_exit (EXIT_FAILURE);
        }
        elem[i].geol.alpha = calib->geol_alpha * geoltbl->alpha[geol_ind];
        elem[i].geol.beta = calib->geol_beta * geoltbl->beta[geol_ind];

        elem[i].geol.dmac = calib->geol_dmac * geoltbl->dmac[geol_ind];
        elem[i].geol.dmac = MIN(elem[i].geol.dmac, elem[i].geol.depth);

        elem[i].geol.areafh = calib->geol_areafh * geoltbl->areafh[geol_ind];
        elem[i].geol.areafv = calib->geol_areafv * geoltbl->areafv[geol_ind];

        elem[i].geol.kmacv = calib->geol_kmacv * geoltbl->kmacv_ro *
            geoltbl->ksatv[geol_ind];
        elem[i].geol.kmach = calib->geol_kmach * geoltbl->kmach_ro *
            geoltbl->ksath[geol_ind];

        elem[i].geol.kinfv  = BADVAL;
        elem[i].geol.dinf   = BADVAL;
        elem[i].geol.smcwlt = BADVAL;
        elem[i].geol.smcref = BADVAL;
#if defined(_NOAH_)
        elem[i].geol.csoil  = BADVAL;
        elem[i].geol.quartz = BADVAL;
        elem[i].geol.smcdry = BADVAL;
#endif
    }
}
