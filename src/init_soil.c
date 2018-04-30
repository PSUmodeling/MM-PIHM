#include "pihm.h"

#if defined(_NOAH_)
void InitSoil(elem_struct *elem, const soiltbl_struct *soiltbl,
    const noahtbl_struct *noahtbl, const calib_struct *cal)
#else
void InitSoil(elem_struct *elem, const soiltbl_struct *soiltbl,
    const calib_struct *cal)
#endif
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
#if defined(_LUMPED_)
    for (i = 0; i < nelem + 1; i++)
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        int             soil_ind;

        soil_ind = elem[i].attrib.soil_type - 1;

        elem[i].soil.dinf = cal->dinf * soiltbl->dinf;

        elem[i].soil.depth = elem[i].topo.zmax - elem[i].topo.zmin;

        elem[i].soil.ksath = cal->ksath * soiltbl->ksath[soil_ind];
        elem[i].soil.ksatv = cal->ksatv * soiltbl->ksatv[soil_ind];
        elem[i].soil.kinfv = cal->kinfv * soiltbl->kinfv[soil_ind];

        elem[i].soil.smcmin = cal->porosity * soiltbl->smcmin[soil_ind];
        elem[i].soil.smcmax = cal->porosity * soiltbl->smcmax[soil_ind];
        elem[i].soil.porosity = elem[i].soil.smcmax - elem[i].soil.smcmin;
        if (elem[i].soil.porosity > 1.0 || elem[i].soil.porosity <= 0.0)
        {
            PIHMprintf(VL_ERROR,
                "Error: Porosity value out of bounds for Element %d", i + 1);
            PIHMexit(EXIT_FAILURE);
        }
        elem[i].soil.alpha = cal->alpha * soiltbl->alpha[soil_ind];
        elem[i].soil.beta = cal->beta * soiltbl->beta[soil_ind];

        /* Calculate field capacity and wilting point following Chan and
         * Dudhia 2001 MWR, but replacing Campbell with van Genuchten */
        elem[i].soil.smcwlt = cal->porosity * soiltbl->smcwlt[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.smcwlt *= cal->smcwlt;
#endif
        elem[i].soil.smcref = cal->porosity * soiltbl->smcref[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.smcref *= cal->smcref;
#endif

        elem[i].soil.dmac = cal->dmac * soiltbl->dmac[soil_ind];
        elem[i].soil.dmac = (elem[i].soil.dmac > elem[i].soil.depth) ?
            elem[i].soil.depth : elem[i].soil.dmac;

        elem[i].soil.areafh = cal->areafh * soiltbl->areafh[soil_ind];
        elem[i].soil.areafv = cal->areafv * soiltbl->areafv[soil_ind];

        elem[i].soil.kmacv =
            cal->kmacv * soiltbl->kmacv_ro * soiltbl->kinfv[soil_ind];
        elem[i].soil.kmach =
            cal->kmach * soiltbl->kmach_ro * soiltbl->ksath[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.csoil = noahtbl->csoil;
        elem[i].soil.quartz = soiltbl->qtz[soil_ind];
        elem[i].soil.smcdry = elem[i].soil.smcwlt;
#endif
    }
}
