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
        int             j, k, kk;
        double          storage;

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;

        for (k = 0; k < NumSpc; k++)
        {
            elem[i].solute[k].conc_surf = elem[i].prcps.t_conc[k] *
                rttbl->Condensation;

            /* Calculate concentrations */
            elem[i].solute[k].conc = (storage > DEPTHR) ?
                elem[i].chms.t_mole[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].solute[k].conc -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms.s_conc[kk];
                    }
                }
            }

            elem[i].solute[k].conc = MAX(elem[i].solute[k].conc, 0.0);
        }

#if defined(_FBR_)
        storage = (elem[i].ws.fbr_unsat + elem[i].ws.fbr_gw) *
            elem[i].geol.porosity + elem[i].geol.depth * elem[i].geol.smcmin;

        for (k = 0; k < NumSpc; k++)
        {
            /* Calculate concentrations */
            elem[i].solute[k].conc_geol = (storage > DEPTHR) ?
                elem[i].chms_geol.t_mole[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].solute[k].conc_geol -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_geol.s_conc[kk];
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
        int             j, k;
        double          storage;

        storage = river[i].ws.stage;

        for (k = 0; k < NumSpc; k++)
        {
            /* Calculate concentrations */
            river[i].solute[k].conc = (storage > DEPTHR) ?
                river[i].chms.t_mole[k] / storage : 0.0;
            river[i].solute[k].conc = MAX(river[i].solute[k].conc, 0.0);
        }
    }
}
