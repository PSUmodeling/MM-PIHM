#include "pihm.h"

#if defined(_NOAH_)
void InitSoil(const soiltbl_struct *soiltbl, const noahtbl_struct *noahtbl, const calib_struct *calib, elem_struct elem[])
#else
void InitSoil(const soiltbl_struct *soiltbl, const calib_struct *calib, elem_struct elem[])
#endif
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             soil_ind;

        if (elem[i].attrib.soil > soiltbl->number)
        {
            pihm_printf(VL_ERROR, "Error: Soil type %d for Element %d is not in the soil file.", elem[i].attrib.soil, i + 1);
            pihm_exit(EXIT_FAILURE);
        }

        soil_ind = elem[i].attrib.soil - 1;

        elem[i].soil.dinf = soiltbl->dinf;

        elem[i].soil.depth = elem[i].topo.zmax - elem[i].topo.zmin;

        elem[i].soil.ksath = calib->ksath * soiltbl->ksath[soil_ind];
        elem[i].soil.ksatv = calib->ksatv * soiltbl->ksatv[soil_ind];
        elem[i].soil.kinfv = calib->kinfv * soiltbl->kinfv[soil_ind];

        elem[i].soil.smcmin = calib->porosity * soiltbl->smcmin[soil_ind];
        elem[i].soil.smcmax = calib->porosity * soiltbl->smcmax[soil_ind];
        elem[i].soil.porosity = elem[i].soil.smcmax - elem[i].soil.smcmin;
        if (elem[i].soil.porosity > 1.0 || elem[i].soil.porosity <= 0.0)
        {
            pihm_printf(VL_ERROR, "Error: Porosity value out of bounds for Element %d", i + 1);
            pihm_exit(EXIT_FAILURE);
        }
        elem[i].soil.alpha = calib->alpha * soiltbl->alpha[soil_ind];
        elem[i].soil.beta = calib->beta * soiltbl->beta[soil_ind];

        // Calculate field capacity and wilting point following Chan and Dudhia 2001 MWR, but replacing Campbell with
        // van Genuchten
        elem[i].soil.smcwlt = calib->porosity * soiltbl->smcwlt[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.smcwlt *= calib->smcwlt;
#endif
        elem[i].soil.smcref = calib->porosity * soiltbl->smcref[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.smcref *= calib->smcref;
#endif

        elem[i].soil.dmac = calib->dmac * soiltbl->dmac[soil_ind];
        elem[i].soil.dmac = MIN(elem[i].soil.dmac, elem[i].soil.depth);

        elem[i].soil.areafh = calib->areafh * soiltbl->areafh[soil_ind];
        elem[i].soil.areafv = calib->areafv * soiltbl->areafv[soil_ind];

        elem[i].soil.kmacv = calib->kmacv * soiltbl->kmacv_ro * soiltbl->kinfv[soil_ind];
        elem[i].soil.kmach = calib->kmach * soiltbl->kmach_ro * soiltbl->ksath[soil_ind];
#if defined(_NOAH_)
        elem[i].soil.csoil = noahtbl->csoil;
        elem[i].soil.quartz = soiltbl->qtz[soil_ind];
        elem[i].soil.smcdry = elem[i].soil.smcwlt;
#endif
    }
}
