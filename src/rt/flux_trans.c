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
            /* Calculate concentrations */
            elem[i].chms.t_conc[k] = (storage > DEPTHR) ?
                elem[i].chms.t_mole[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].chms.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms.s_conc[kk];
                    }
                }
            }

            elem[i].chms.t_conc[k] = MAX(elem[i].chms.t_conc[k], 0.0);
        }

#if defined(_FBR_)
        storage = (elem[i].ws.fbr_unsat + elem[i].ws.fbr_gw) *
            elem[i].geol.porosity + elem[i].geol.depth * elem[i].geol.smcmin;

        for (k = 0; k < NumSpc; k++)
        {
            /* Calculate concentrations */
            elem[i].chms_geol.t_conc[k] = (storage > DEPTHR) ?
                elem[i].chms_geol.t_mole[k] / storage : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].chms_geol.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_geol.s_conc[kk];
                    }
                }
            }

            elem[i].chms_geol.t_conc[k] = MAX(elem[i].chms_geol.t_conc[k], 0.0);
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
            river[i].chms.t_conc[k] = (storage > DEPTHR) ?
                river[i].chms.t_mole[k] / storage : 0.0;
            river[i].chms.t_conc[k] = MAX(river[i].chms.t_conc[k], 0.0);
        }
    }
}
