#include "pihm.h"

void InitLc(const lctbl_struct *lctbl, const calib_struct *calib, elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        _InitLc(lctbl, calib, &elem[i]);
    }
}

void _InitLc(const lctbl_struct *lctbl, const calib_struct *calib, elem_struct *elem_ptr)
{
    int             lc_ind;

    if (elem_ptr->attrib.lc > lctbl->number)
    {
        pihm_printf(VL_ERROR, "Error: Land cover type %d for Element %d is not in the vegetation file.",
            elem_ptr->attrib.lc, elem_ptr->ind);
        pihm_exit(EXIT_FAILURE);
    }

    lc_ind = elem_ptr->attrib.lc - 1;

    elem_ptr->ps.rzd = calib->rzd * lctbl->rzd[lc_ind];
    elem_ptr->epc.rsmin = lctbl->rsmin[lc_ind];
    elem_ptr->epc.rgl = lctbl->rgl[lc_ind];
    elem_ptr->epc.hs = lctbl->hs[lc_ind];
    elem_ptr->epc.rsmax = lctbl->rsmax;
    elem_ptr->epc.topt = lctbl->topt;
    elem_ptr->lc.shdfac = calib->vegfrac * lctbl->vegfrac[lc_ind];
    elem_ptr->lc.laimin = lctbl->laimin[lc_ind];
    elem_ptr->lc.laimax = lctbl->laimax[lc_ind];
    elem_ptr->lc.emissmin = lctbl->emissmin[lc_ind];
    elem_ptr->lc.emissmax = lctbl->emissmax[lc_ind];
    elem_ptr->lc.albedomin = calib->albedo * lctbl->albedomin[lc_ind];
    elem_ptr->lc.albedomax = calib->albedo * lctbl->albedomax[lc_ind];
    elem_ptr->lc.z0min = lctbl->z0min[lc_ind];
    elem_ptr->lc.z0max = lctbl->z0max[lc_ind];
    elem_ptr->lc.rough = calib->rough * lctbl->rough[lc_ind];
    elem_ptr->lc.cmcfactr = CMCFACTR;
    elem_ptr->lc.cfactr = lctbl->cfactr;
    elem_ptr->lc.bare = (IGBP_BARREN == elem_ptr->attrib.lc || NLCD40_BARREN == elem_ptr->attrib.lc) ? 1 : 0;
    elem_ptr->lc.shdfac = (1 == elem_ptr->lc.bare) ? 0.0 : elem_ptr->lc.shdfac;
#if defined(_NOAH_)
    elem_ptr->lc.snup = lctbl->snup[lc_ind];
    elem_ptr->lc.isurban = (IGBP_URBAN_BUILDUP == elem_ptr->attrib.lc || NLCD40_DEVELOPED_OPEN == elem_ptr->attrib.lc || NLCD40_DEVELOPED_LOW == elem_ptr->attrib.lc ||
        NLCD40_DEVELOPED_MID == elem_ptr->attrib.lc || NLCD40_DEVELOPED_HIGH == elem_ptr->attrib.lc) ? 1 : 0;
    elem_ptr->lc.glacier = (IGBP_SNOW_ICE == elem_ptr->attrib.lc || NLCD40_SNOW_ICE == elem_ptr->attrib.lc) ? 1 : 0;
    elem_ptr->lc.shdmin = 0.01;
    elem_ptr->lc.shdmax = 0.96;
#endif

#if defined(_NOAH_)
    elem_ptr->epc.rgl *= calib->rgl;
    elem_ptr->epc.hs *= calib->hs;
    elem_ptr->epc.rsmin *= calib->rsmin;
    elem_ptr->lc.cmcfactr *= calib->cmcmax;
    elem_ptr->lc.cfactr *= calib->cfactr;
#endif
}
