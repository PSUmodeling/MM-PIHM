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
        int             lc_ind;

        lc_ind = elem[i].attrib.lc_type - 1;

        elem[i].ps.rzd = cal->rzd * lctbl->rzd[lc_ind];
        elem[i].epc.rsmin = lctbl->rsmin[lc_ind];
        elem[i].epc.rgl = lctbl->rgl[lc_ind];
        elem[i].epc.hs = lctbl->hs[lc_ind];
        elem[i].epc.rsmax = lctbl->rsmax;
        elem[i].epc.topt = lctbl->topt;
        elem[i].lc.shdfac = cal->vegfrac * lctbl->vegfrac[lc_ind];
        elem[i].lc.laimin = lctbl->laimin[lc_ind];
        elem[i].lc.laimax = lctbl->laimax[lc_ind];
        elem[i].lc.emissmin = lctbl->emissmin[lc_ind];
        elem[i].lc.emissmax = lctbl->emissmax[lc_ind];
        elem[i].lc.albedomin = cal->albedo * lctbl->albedomin[lc_ind];
        elem[i].lc.albedomax = cal->albedo * lctbl->albedomax[lc_ind];
        elem[i].lc.z0min = lctbl->z0min[lc_ind];
        elem[i].lc.z0max = lctbl->z0max[lc_ind];
        elem[i].lc.rough = cal->rough * lctbl->rough[lc_ind];
        elem[i].lc.cmcfactr = CMCFACTR;
        elem[i].lc.cfactr = lctbl->cfactr;
        elem[i].lc.bare = (lctbl->bare == elem[i].attrib.lc_type) ? 1 : 0;
        elem[i].lc.shdfac = (1 == elem[i].lc.bare) ? 0.0 : elem[i].lc.shdfac;
#if defined(_NOAH_)
        elem[i].lc.snup = lctbl->snup[lc_ind];
        elem[i].lc.isurban = (ISURBAN == elem[i].attrib.lc_type) ? 1 : 0;
        elem[i].lc.shdmin = 0.01;
        elem[i].lc.shdmax = 0.96;
#endif

#if defined(_NOAH_)
        elem[i].epc.rgl *= cal->rgl;
        elem[i].epc.hs *= cal->hs;
        elem[i].epc.rsmin *= cal->rsmin;
        elem[i].lc.cmcfactr *= cal->cmcmax;
        elem[i].lc.cfactr *= cal->cfactr;
#endif
    }
}
