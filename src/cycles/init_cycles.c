#include "pihm.h"

void InitCycles(const agtbl_struct *agtbl, const mgmt_struct mgmttbl[],
    const crop_struct croptbl[], const soiltbl_struct *soiltbl,
    elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             soil_ind;
        int             k, kcrop;
        /*
         * Initialize initial soil variables
         */
        soil_ind = elem[i].attrib.soil_type - 1;

        for (k = 0; k < MAXLYR; k++)
        {
            if (k < elem[i].ps.nsoil)
            {
                elem[i].soil.clay[k]  = soiltbl->clay_layer[soil_ind][k];
                elem[i].soil.sand[k]  = soiltbl->sand_layer[soil_ind][k];
                elem[i].soil.om[k]    = soiltbl->om_layer[soil_ind][k];
                elem[i].soil.bd[k]    = soiltbl->bd_layer[soil_ind][k];
                elem[i].soil.fc[k]    = soiltbl->fc[soil_ind][k];
                elem[i].soil.pwp[k]   = soiltbl->pwp[soil_ind][k];
                elem[i].soil.b[k]     = soiltbl->b[soil_ind][k];
                elem[i].soil.air_entry_pot[k] =
                                        soiltbl->air_entry_pot[soil_ind][k];
            }
            else
            {
                elem[i].soil.clay[k]          = BADVAL;
                elem[i].soil.sand[k]          = BADVAL;
                elem[i].soil.om[k]            = BADVAL;
                elem[i].soil.bd[k]            = BADVAL;
                elem[i].soil.fc[k]            = BADVAL;
                elem[i].soil.pwp[k]           = BADVAL;
                elem[i].soil.b[k]             = BADVAL;
                elem[i].soil.air_entry_pot[k] = BADVAL;
            }
        }

        /*
         * Initialize crop structures
         */
        for (kcrop = 0; kcrop < MAXCROP; kcrop++)
        {
            elem[i].crop[kcrop] = croptbl[kcrop];

            if (croptbl[kcrop].epc.name[0] != '\0')
            {
                InitCropStateVar(&elem[i].crop[kcrop]);
            }
        }

        /*
         * Initialize management structure
         */
        elem[i].mgmt = mgmttbl[agtbl->oper[i] - 1];
    }
}

#if defined(_CYCLES_OBSOLETE_)
void FirstDay(const soiltbl_struct *soiltbl, elem_struct elem[],
    river_struct river[])
{
    int             i, k;
    int             soil_ind;

    for (i = 0; i < nelem; i++)
    {
        soil_ind = elem[i].attrib.soil_type - 1;

        elem[i].restart_input.resw_stan = 0.0;
        elem[i].restart_input.resw_flat = 0.0;
        elem[i].restart_input.resm_stan = 0.0;
        elem[i].restart_input.resm_flat = 0.0;
        elem[i].restart_input.manuc_surf = 0.0;
        elem[i].restart_input.resn_stan = 0.0;
        elem[i].restart_input.resn_flat = 0.0;
        elem[i].restart_input.manun_surf = 0.0;

        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].restart_input.res_abgd[k] = BADVAL;
            elem[i].restart_input.res_root[k] = BADVAL;
            elem[i].restart_input.res_rhizo[k] = BADVAL;
            elem[i].restart_input.manuc[k] = BADVAL;
            elem[i].restart_input.resn_abgd[k] = BADVAL;
            elem[i].restart_input.resn_root[k] = BADVAL;
            elem[i].restart_input.resn_rhizo[k] = BADVAL;
            elem[i].restart_input.manun[k] = BADVAL;
            elem[i].restart_input.soc[k] = BADVAL;
            elem[i].restart_input.son[k] = BADVAL;
            elem[i].restart_input.mbc[k] = BADVAL;
            elem[i].restart_input.mbn[k] = BADVAL;
            elem[i].restart_input.no3[k] = BADVAL;
            elem[i].restart_input.nh4[k] = BADVAL;
        }

        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            elem[i].restart_input.soc[k] = elem[i].soil.iom[k] / 100.0 * 0.58 *
                elem[i].ps.sldpth[k] * elem[i].soil.bd[k] * 1000.0;
            /* Initializes as 3% of SOC_Mass but "added" C */
            elem[i].restart_input.mbc[k] = 0.03 * elem[i].restart_input.soc[k];
            elem[i].restart_input.res_abgd[k] = 0.0;
            elem[i].restart_input.res_root[k] = 0.0;
            elem[i].restart_input.res_rhizo[k] = 0.0;
            elem[i].restart_input.manuc[k] = 0.0;

            elem[i].restart_input.no3[k] =
                soiltbl->no3_lyr[soil_ind][k] * 1.0E-4;
            elem[i].restart_input.nh4[k] =
                soiltbl->nh4_lyr[soil_ind][k] * 1.0E-4;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.son[k] = elem[i].restart_input.soc[k] * 0.1;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.mbn[k] = elem[i].restart_input.mbc[k] * 0.1;
            elem[i].restart_input.resn_abgd[k] = 0.0;
            elem[i].restart_input.resn_root[k] = 0.0;
            elem[i].restart_input.resn_rhizo[k] = 0.0;
            elem[i].restart_input.manun[k] = 0.0;
        }
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].restart_input.streamno3 = 0.0;
        river[i].restart_input.streamnh4 = 0.0;
        river[i].restart_input.bedno3 = 0.0;
        river[i].restart_input.bednh4 = 0.0;
    }
}

void InitCyclesVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i, k;

    for (i = 0; i < nelem; i++)
    {
        RestartInput(&elem[i].restart_input, &elem[i].ps, &elem[i].ws,
            &elem[i].cs, &elem[i].ns);

        MakeZeroFluxStruct(&elem[i].wf, &elem[i].cf, &elem[i].nf);
        elem[i].no3sol.snksrc = 0.0;
        elem[i].nh4sol.snksrc = 0.0;

        elem[i].np.no3 = 0.0;
        elem[i].np.nh4 = 0.0;
        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            elem[i].np.no3 += elem[i].ns.no3[k];
            elem[i].np.nh4 += elem[i].ns.nh4[k];
        }

        NV_Ith(CV_Y, NO3(i)) = elem[i].np.no3;
        NV_Ith(CV_Y, NH4(i)) = elem[i].np.nh4;

        elem[i].ns0 = elem[i].ns;
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].ns.streamno3 = river[i].restart_input.streamno3;
        river[i].ns.streamnh4 = river[i].restart_input.streamnh4;
        river[i].ns.bedno3 = river[i].restart_input.bedno3;
        river[i].ns.bednh4 = river[i].restart_input.bednh4;

        NV_Ith(CV_Y, STREAMNO3(i)) = river[i].ns.streamno3;
        NV_Ith(CV_Y, STREAMNH4(i)) = river[i].ns.streamnh4;
        NV_Ith(CV_Y, RIVBEDNO3(i)) = river[i].ns.bedno3;
        NV_Ith(CV_Y, RIVBEDNH4(i)) = river[i].ns.bednh4;
    }
}
//
//void WriteCyclesIC(char *restart_fn, elem_struct *elem, river_struct *riv)
//{
//    int             i, j;
//    FILE           *restart_file;
//
//    restart_file = PIHMfopen(restart_fn, "wb");
//    PIHMprintf(VL_VERBOSE, "Writing Cycles initial conditions.\n");
//
//    for (i = 0; i < nelem; i++)
//    {
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.SOC_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.SON_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.MBC_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.MBN_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.NO3[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.NH4[j], sizeof(double), 1, restart_file);
//        }
//    }
//
//    for (i = 0; i < nriver; i++)
//    {
//        fwrite(&riv[i].NO3sol.soluteMass, sizeof(double), 1, restart_file);
//        fwrite(&riv[i].NH4sol.soluteMass, sizeof(double), 1, restart_file);
//    }
//
//    fclose(restart_file);
//}
#endif