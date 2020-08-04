#include "pihm.h"

void SoluteConc(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[])
{
    int             i;

    /*
     * Calculate chemical concentrations
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k, kk;
        double          storage;

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].conc_surf = elem[i].prcpchm.tot_conc[k] *
                rttbl->cond;

            /* Calculate concentrations */
            elem[i].solute[k].conc = (storage > DEPTHR) ?
                elem[i].chms.tot_mol[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->num_ssc; kk++)
                {
                    if ((rttbl->conc_contrib[k][kk + rttbl->num_stc] != 0) &&
                        (chemtbl[kk + rttbl->num_stc].itype != AQUEOUS))
                    {
                        elem[i].solute[k].conc -=
                            rttbl->conc_contrib[k][kk + rttbl->num_stc] *
                            elem[i].chms.sec_conc[kk];
                    }
                }
            }

            elem[i].solute[k].conc = MAX(elem[i].solute[k].conc, 0.0);
        }

#if defined(_DGW_)
        storage = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) *
            elem[i].geol.porosity + elem[i].geol.depth * elem[i].geol.smcmin;

        for (k = 0; k < nsolute; k++)
        {
            /* Calculate concentrations */
            elem[i].solute[k].conc_geol = (storage > DEPTHR) ?
                elem[i].chms_geol.tot_mol[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->num_ssc; kk++)
                {
                    if ((rttbl->conc_contrib[k][kk + rttbl->num_stc] != 0) &&
                        (chemtbl[kk + rttbl->num_stc].itype != AQUEOUS))
                    {
                        elem[i].solute[k].conc_geol -=
                            rttbl->conc_contrib[k][kk + rttbl->num_stc] *
                            elem[i].chms_geol.sec_conc[kk];
                    }
                }
            }

            elem[i].solute[k].conc_geol = MAX(elem[i].solute[k].conc_geol, 0.0);
        }
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;
        double          storage;

        storage = river[i].ws.stage;

        for (k = 0; k < nsolute; k++)
        {
            /* Calculate concentrations */
            river[i].solute[k].conc = (storage > DEPTHR) ?
                river[i].chms.tot_mol[k] / storage : 0.0;
            river[i].solute[k].conc = MAX(river[i].solute[k].conc, 0.0);
        }
    }
}
