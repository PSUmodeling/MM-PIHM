#include "pihm.h"

void InitLc(elem_struct *elem, const lctbl_struct *lctbl,
    const calib_struct *cal)
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
        _InitLc(&elem[i], lctbl, cal);
    }
}

void _InitLc(elem_struct *elem_ptr, const lctbl_struct *lctbl,
    const calib_struct *cal)
{
    int             lc_ind;

    lc_ind = elem_ptr->attrib.lc_type - 1;

    elem_ptr->ps.rzd = cal->rzd * lctbl->rzd[lc_ind];
#if !defined(_CYCLES_)
    elem_ptr->epc.rsmin = lctbl->rsmin[lc_ind];
    elem_ptr->epc.rgl = lctbl->rgl[lc_ind];
    elem_ptr->epc.hs = lctbl->hs[lc_ind];
    elem_ptr->epc.rsmax = lctbl->rsmax;
    elem_ptr->epc.topt = lctbl->topt;
#endif
    elem_ptr->lc.shdfac = cal->vegfrac * lctbl->vegfrac[lc_ind];
    elem_ptr->lc.laimin = lctbl->laimin[lc_ind];
    elem_ptr->lc.laimax = lctbl->laimax[lc_ind];
    elem_ptr->lc.emissmin = lctbl->emissmin[lc_ind];
    elem_ptr->lc.emissmax = lctbl->emissmax[lc_ind];
    elem_ptr->lc.albedomin = cal->albedo * lctbl->albedomin[lc_ind];
    elem_ptr->lc.albedomax = cal->albedo * lctbl->albedomax[lc_ind];
    elem_ptr->lc.z0min = lctbl->z0min[lc_ind];
    elem_ptr->lc.z0max = lctbl->z0max[lc_ind];
    elem_ptr->lc.rough = cal->rough * lctbl->rough[lc_ind];
    elem_ptr->lc.cmcfactr = CMCFACTR;
    elem_ptr->lc.cfactr = lctbl->cfactr;
    elem_ptr->lc.bare =
        (IGBP_BARREN == elem_ptr->attrib.lc_type ||
        NLCD40_BARREN == elem_ptr->attrib.lc_type) ? 1 : 0;
    elem_ptr->lc.shdfac = (1 == elem_ptr->lc.bare) ? 0.0 : elem_ptr->lc.shdfac;
#if defined(_NOAH_)
    elem_ptr->lc.snup = lctbl->snup[lc_ind];
    elem_ptr->lc.isurban =
        (IGBP_URBAN_BUILDUP == elem_ptr->attrib.lc_type ||
        NLCD40_DEVELOPED_OPEN == elem_ptr->attrib.lc_type ||
        NLCD40_DEVELOPED_LOW == elem_ptr->attrib.lc_type ||
        NLCD40_DEVELOPED_MID == elem_ptr->attrib.lc_type ||
        NLCD40_DEVELOPED_HIGH == elem_ptr->attrib.lc_type) ? 1 : 0;
    elem_ptr->lc.glacier =
        (IGBP_SNOW_ICE == elem_ptr->attrib.lc_type ||
        NLCD40_SNOW_ICE == elem_ptr->attrib.lc_type) ? 1 : 0;
    elem_ptr->lc.shdmin = 0.01;
    elem_ptr->lc.shdmax = 0.96;
#endif

#if defined(_NOAH_)
# if !defined(_CYCLES_)
    elem_ptr->epc.rgl *= cal->rgl;
    elem_ptr->epc.hs *= cal->hs;
    elem_ptr->epc.rsmin *= cal->rsmin;
# endif
    elem_ptr->lc.cmcfactr *= cal->cmcmax;
    elem_ptr->lc.cfactr *= cal->cfactr;
#endif
}
