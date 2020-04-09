/*******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that operates
* on top of the hydrological land surface processes described by Flux-PIHM.
* RT-Flux-PIHM tracks the transportation and reaction in a given watershed. It
* uses operator splitting technique to couple transport and reaction.
*****************************************************************************/
#include "pihm.h"

void InitChem(const char cdbs_filen[], const calib_struct *cal,
    forc_struct *forc, chemtbl_struct chemtbl[], kintbl_struct kintbl[],
    rttbl_struct *rttbl, chmictbl_struct *chmictbl, elem_struct elem[])
{
    int             i, j;
    int             chem_ind;
    FILE           *fp;

    fp = fopen(cdbs_filen, "r");
    CheckFile(fp, cdbs_filen);

    /*
     * Look up database to find required parameters and dependencies for
     * chemical species
     */
    Lookup(fp, cal, chemtbl, kintbl, rttbl);
    fclose(fp);

    /*
     * Apply calibration
     */
    chem_ind = FindChem("'DOC'", chemtbl, rttbl->NumStc);
    if (chem_ind >= 0)
    {
        rttbl->prcp_conc[chem_ind] *= cal->prcpconc;

        if (forc->PrpFlg == 2)
        {
            for (i = 0; i < forc->nprcpc; i++)
            {
                for (j = 0; j < forc->prcpc[i].length; j++)
                {
                    forc->prcpc[i].data[j][chem_ind] *= cal->prcpconc;
                }
            }
        }
    }

    for (i = 0; i < chmictbl->nic; i++)
    {
        int             k;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            chmictbl->ssa[i][k] *= (chemtbl[k].itype == MINERAL) ?
                cal->ssa : 1.0;

            chmictbl->conc[i][k] *=
                (strcmp(chemtbl[k].ChemName, "'DOC'") == 0) ?
                cal->initconc : 1.0;
        }
    }

    /*
     * Assign initial conditions to different volumes
     */
    for (i = 0; i < nelem; i++)
    {
        int             k;
        int            *ic_type;

        ic_type = elem[i].attrib.chem_ic_type;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            elem[i].restart_input[SOIL_CHMVOL].t_conc[k] =
                chmictbl->conc[ic_type[SOIL_CHMVOL] - 1][k];
            elem[i].restart_input[SOIL_CHMVOL].ssa[k] =
                chmictbl->ssa[ic_type[SOIL_CHMVOL] - 1][k];

#if defined(_FBR_)
            elem[i].restart_input[GEOL_CHMVOL].t_conc[k] =
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
    PIHMprintf(VL_VERBOSE, "\n Initializing concentrations... \n");

    for (i = 0; i < nelem; i++)
    {
        double          storage;

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;

        InitChemS(chemtbl, rttbl, &elem[i].restart_input[SOIL_CHMVOL],
            elem[i].soil.smcmax, storage, &elem[i].chms);

#if defined(_FBR_)
        storage = (elem[i].ws.fbr_unsat + elem[i].ws.fbr_gw) *
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

        for (k = 0; k < rttbl->NumStc; k++)
        {
            if (chemtbl[k].itype == AQUEOUS)
            {
                river[i].chms.t_conc[k] =
                    0.5 * elem[river[i].leftele - 1].chms.t_conc[k] +
                    0.5 * elem[river[i].rightele - 1].chms.t_conc[k];
                river[i].chms.p_actv[k] = river[i].chms.t_conc[k];
                river[i].chms.p_conc[k] = river[i].chms.t_conc[k];
                river[i].chms.t_mole[k] = river[i].chms.t_conc[k] * storage;
            }
            else
            {
                river[i].chms.t_conc[k] = ZERO_CONC;
                river[i].chms.p_conc[k] = ZERO_CONC;
                river[i].chms.p_actv[k] = ZERO_CONC;
                river[i].chms.t_mole[k] = 0.0;
            }
        }

        for (k = 0; k < rttbl->NumSsc; k++)
        {
            river[i].chms.s_conc[k] = ZERO_CONC;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, SOIL_MOLE(i, k)) = elem[i].chms.t_mole[k];

            elem[i].chmf.react[k] = 0.0;

#if defined(_FBR_)
            NV_Ith(CV_Y, GEOL_MOLE(i, k)) = elem[i].chms_geol.t_mole[k];

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

        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, RIVER_MOLE(i, k)) = river[i].chms.t_mole[k];
        }
    }
}

void InitChemS(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    const rtic_struct *restart_input, double smcmax, double vol,
    chmstate_struct *chms)
{
    int             k;

    for (k = 0; k < rttbl->NumStc; k++)
    {
        if (strcmp(chemtbl[k].ChemName, "'H+'") == 0)
        {
            chms->t_conc[k] = restart_input->t_conc[k];
            chms->p_actv[k] = chms->t_conc[k];
            chms->p_conc[k] = chms->t_conc[k];
            chms->ssa[k] = restart_input->ssa[k];
        }
        else if (chemtbl[k].itype == MINERAL)
        {
            chms->t_conc[k] = restart_input->t_conc[k];
            /* Update the concentration of mineral using molar volume */
            chms->t_conc[k] *= (rttbl->RelMin == 0) ?
                /* Absolute mineral volume fraction */
                1000.0 / chemtbl[k].MolarVolume / smcmax :
                /* Relative mineral volume fraction */
                (1.0 - smcmax) * 1000.0 / chemtbl[k].MolarVolume / smcmax;
            chms->p_actv[k] = 1.0;
            chms->p_conc[k] = chms->t_conc[k];
            chms->ssa[k] = restart_input->ssa[k];
        }
        else if ((chemtbl[k].itype == CATION_ECHG) ||
            (chemtbl[k].itype == ADSORPTION))
        {
            chms->t_conc[k] = restart_input->t_conc[k];
            chms->p_actv[k] = chms->t_conc[k] * 0.5;
            /* Change unit of CEC (eq g-1) into C(ion site)
             * (eq L-1 porous space), assuming density of solid is always
             * 2650 g L-1 */
            chms->t_conc[k] *= (1.0 - smcmax) * 2650.0;
            chms->p_conc[k] = chms->t_conc[k];
        }
        else
        {
            chms->t_conc[k] = restart_input->t_conc[k];
            chms->p_actv[k] = chms->t_conc[k] * 0.5;
            chms->p_conc[k] = chms->t_conc[k] * 0.5;
            chms->ssa[k] = restart_input->ssa[k];
        }
    }

    for (k = 0; k < rttbl->NumSsc; k++)
    {
        chms->s_conc[k] = ZERO_CONC;
    }

    /* Speciation */
    if (rttbl->RecFlg == KIN_REACTION)
    {
        _Speciation(chemtbl, rttbl, 1, chms);
    }

    /* Total moles should be calculated after speciation */
    for (k = 0; k < rttbl->NumStc; k++)
    {
        if (chemtbl[k].itype == AQUEOUS)
        {
            chms->t_mole[k] = chms->t_conc[k] * vol;
        }
        else
        {
            chms->t_mole[k] = 0.0;
        }
    }
}

void UpdatePConc(elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            elem[i].chms.p_conc[k] = elem[i].chms.t_conc[k];

#if defined(_FBR_)
            elem[i].chms_geol.p_conc[k] = elem[i].chms_geol.t_conc[k];
#endif
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            river[i].chms.p_conc[k] = river[i].chms.t_conc[k];
        }
    }
}
