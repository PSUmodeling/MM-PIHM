/*******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that operates
* on top of the hydrological land surface processes described by Flux-PIHM.
* RT-Flux-PIHM tracks the transportation and reaction in a given watershed. It
* uses operator splitting technique to couple transport and reaction.
*****************************************************************************/
#include "pihm.h"

void InitChem(const char cdbs_filen[], const calib_struct *calib,
    forc_struct *forc, chemtbl_struct chemtbl[], kintbl_struct kintbl[],
    rttbl_struct *rttbl, chmictbl_struct *chmictbl, elem_struct elem[])
{
    int             i, j;
    int             chem_ind;
    FILE           *fp;

    fp = pihm_fopen(cdbs_filen, "r");

    /*
     * Look up database to find required parameters and dependencies for
     * chemical species
     */
    Lookup(fp, calib, chemtbl, kintbl, rttbl);
    fclose(fp);

    /*
     * Apply calibration
     */
    chem_ind = FindChem("'DOC'", chemtbl, rttbl->num_stc);
    if (chem_ind >= 0)
    {
        rttbl->prcp_conc[chem_ind] *= calib->prcpconc;

        if (forc->prcp_flag == 2)
        {
            for (i = 0; i < forc->nprcpc; i++)
            {
                for (j = 0; j < forc->prcpc[i].length; j++)
                {
                    forc->prcpc[i].data[j][chem_ind] *= calib->prcpconc;
                }
            }
        }
    }

    for (i = 0; i < chmictbl->nic; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_stc; k++)
        {
            chmictbl->ssa[i][k] *= (chemtbl[k].itype == MINERAL) ?
                calib->ssa : 1.0;

            chmictbl->conc[i][k] *=
                (strcmp(chemtbl[k].name, "'DOC'") == 0) ?
                calib->initconc : 1.0;
        }
    }

    /*
     * Assign initial conditions to different volumes
     */
    for (i = 0; i < nelem; i++)
    {
        int             k;
        int            *ic_type;

        ic_type = elem[i].attrib.chem_ic;

        for (k = 0; k < rttbl->num_stc; k++)
        {
            elem[i].restart_input[SOIL_CHMVOL].tot_conc[k] =
                chmictbl->conc[ic_type[SOIL_CHMVOL] - 1][k];
            elem[i].restart_input[SOIL_CHMVOL].ssa[k] =
                chmictbl->ssa[ic_type[SOIL_CHMVOL] - 1][k];

#if defined(_DGW_)
            elem[i].restart_input[GEOL_CHMVOL].tot_conc[k] =
                chmictbl->conc[ic_type[GEOL_CHMVOL] - 1][k];
            elem[i].restart_input[GEOL_CHMVOL].ssa[k] =
                chmictbl->ssa[ic_type[GEOL_CHMVOL] - 1][k];
#endif
        }
    }
}

void InitRTVar(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

    /*
     * Initializing element concentrations
     */
    pihm_printf(VL_VERBOSE, "\n Initializing concentrations... \n");

    for (i = 0; i < nelem; i++)
    {
        double          storage;

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;

        InitChemS(chemtbl, rttbl, &elem[i].restart_input[SOIL_CHMVOL],
            elem[i].soil.smcmax, storage, &elem[i].chms);

#if defined(_DGW_)
        storage = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) *
            elem[i].geol.porosity + elem[i].geol.depth * elem[i].geol.smcmin;

        InitChemS(chemtbl, rttbl, &elem[i].restart_input[GEOL_CHMVOL],
            elem[i].geol.smcmax, storage, &elem[i].chms_geol);
#endif
    }

    /*
     * Initialize river concentrations
     */
    for (i = 0; i < nriver; i++)
    {
        double          storage;
        int             k;

        storage = MAX(river[i].ws.stage, DEPTHR);

        for (k = 0; k < rttbl->num_stc; k++)
        {
            if (chemtbl[k].itype == AQUEOUS)
            {
                river[i].chms.tot_conc[k] =
                    0.5 * elem[river[i].left - 1].chms.tot_conc[k] +
                    0.5 * elem[river[i].right - 1].chms.tot_conc[k];
                river[i].chms.prim_actv[k] = river[i].chms.tot_conc[k];
                river[i].chms.prim_conc[k] = river[i].chms.tot_conc[k];
                river[i].chms.tot_mol[k] = river[i].chms.tot_conc[k] * storage;
            }
            else
            {
                river[i].chms.tot_conc[k] = ZERO_CONC;
                river[i].chms.prim_conc[k] = ZERO_CONC;
                river[i].chms.prim_actv[k] = ZERO_CONC;
                river[i].chms.tot_mol[k] = 0.0;
            }
        }

        for (k = 0; k < rttbl->num_ssc; k++)
        {
            river[i].chms.sec_conc[k] = ZERO_CONC;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_spc; k++)
        {
            NV_Ith(CV_Y, SOLUTE_SOIL(i, k)) = elem[i].chms.tot_mol[k];

            elem[i].chmf.react[k] = 0.0;

#if defined(_DGW_)
            NV_Ith(CV_Y, SOLUTE_GEOL(i, k)) = elem[i].chms_geol.tot_mol[k];

            elem[i].chmf.react_geol[k] = 0.0;
#endif
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_spc; k++)
        {
            NV_Ith(CV_Y, SOLUTE_RIVER(i, k)) = river[i].chms.tot_mol[k];
        }
    }
}

void InitChemS(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    const rtic_struct *restart_input, double smcmax, double vol,
    chmstate_struct *chms)
{
    int             k;

    for (k = 0; k < rttbl->num_stc; k++)
    {
        if (strcmp(chemtbl[k].name, "'H+'") == 0)
        {
            chms->tot_conc[k] = restart_input->tot_conc[k];
            chms->prim_actv[k] = chms->tot_conc[k];
            chms->prim_conc[k] = chms->tot_conc[k];
            chms->ssa[k] = restart_input->ssa[k];
        }
        else if (chemtbl[k].itype == MINERAL)
        {
            chms->tot_conc[k] = restart_input->tot_conc[k];
            /* Update the concentration of mineral using molar volume */
            chms->tot_conc[k] *= (rttbl->rel_min == 0) ?
                /* Absolute mineral volume fraction */
                1000.0 / chemtbl[k].molar_vol / smcmax :
                /* Relative mineral volume fraction */
                (1.0 - smcmax) * 1000.0 / chemtbl[k].molar_vol / smcmax;
            chms->prim_actv[k] = 1.0;
            chms->prim_conc[k] = chms->tot_conc[k];
            chms->ssa[k] = restart_input->ssa[k];
        }
        else if ((chemtbl[k].itype == CATION_ECHG) ||
            (chemtbl[k].itype == ADSORPTION))
        {
            chms->tot_conc[k] = restart_input->tot_conc[k];
            chms->prim_actv[k] = chms->tot_conc[k] * 0.5;
            /* Change unit of CEC (eq g-1) into C(ion site)
             * (eq L-1 porous space), assuming density of solid is always
             * 2650 g L-1 */
            chms->tot_conc[k] *= (1.0 - smcmax) * 2650.0;
            chms->prim_conc[k] = chms->tot_conc[k];
        }
        else
        {
            chms->tot_conc[k] = restart_input->tot_conc[k];
            chms->prim_actv[k] = chms->tot_conc[k] * 0.5;
            chms->prim_conc[k] = chms->tot_conc[k] * 0.5;
            chms->ssa[k] = restart_input->ssa[k];
        }
    }

    for (k = 0; k < rttbl->num_ssc; k++)
    {
        chms->sec_conc[k] = ZERO_CONC;
    }

    /* Speciation */
    if (rttbl->transpt_flag == KIN_REACTION)
    {
        _Speciation(chemtbl, rttbl, 1, chms);
    }

    /* Total moles should be calculated after speciation */
    for (k = 0; k < rttbl->num_stc; k++)
    {
        if (chemtbl[k].itype == AQUEOUS)
        {
            chms->tot_mol[k] = chms->tot_conc[k] * vol;
        }
        else
        {
            chms->tot_mol[k] = 0.0;
        }
    }
}

void UpdatePConc(const rttbl_struct *rttbl, elem_struct elem[],
    river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_spc; k++)
        {
            elem[i].chms.prim_conc[k] = elem[i].chms.tot_conc[k];

#if defined(_DGW_)
            elem[i].chms_geol.prim_conc[k] = elem[i].chms_geol.tot_conc[k];
#endif
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_spc; k++)
        {
            river[i].chms.prim_conc[k] = river[i].chms.tot_conc[k];
        }
    }
}
