#include "pihm.h"

void InitCycles(const agtbl_struct *agtbl, const soiltbl_struct *soiltbl,
    epconst_struct epctbl[], elem_struct elem[], river_struct river[])
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        int             soil_ind;
        int             j, k;
        /*
         * Initialize initial soil variables
         */
        soil_ind = elem[i].attrib.soil_type - 1;

        for (k = 0; k < MAXLYR; k++)
        {
            if (k < elem[i].ps.nsoil)
            {
                elem[i].soil.clay[k] = soiltbl->clay_lyr[soil_ind][k];
                elem[i].soil.clay[k] *= 0.01;
                elem[i].soil.sand[k] = soiltbl->sand_lyr[soil_ind][k];
                elem[i].soil.sand[k] *= 0.01;
                elem[i].soil.iom[k] = soiltbl->iom_lyr[soil_ind][k];
                elem[i].soil.bd[k] = soiltbl->bd_lyr[soil_ind][k];
            }
            else
            {
                elem[i].soil.clay[k] = BADVAL;
                elem[i].soil.sand[k] = BADVAL;
                elem[i].soil.iom[k] = BADVAL;
                elem[i].soil.bd[k] = BADVAL;
            }
        }

        /*
         * Initialize crop structures
         */
        for (j = 0; j < MAXCROP; j++)
        {
            elem[i].crop[j].epc = &epctbl[j];

            InitCropSV(&elem[i].crop[j]);
        }

        /*
         * Initialize management structure
         */
        elem[i].mgmt.rot_size = agtbl->rotsz[i];
        elem[i].mgmt.auto_n = agtbl->auto_N[i];
        elem[i].mgmt.rot_year = 0;

        for (j = 0; j < 4; j++)
        {
            elem[i].mgmt.op_ptr[j] = 0;
        }
    }

    /*
     * Initialize river bed bulk densities
     */
    for (i = 0; i < nriver; i++)
    {
        int             k;
        double          totdpth;

        river[i].matl.bd = 0.0;
        totdpth = 0.0;

        for (k = 0; k < elem[river[i].leftele - 1].ps.nsoil; k++)
        {
            river[i].matl.bd += elem[river[i].leftele - 1].soil.bd[k] *
                elem[river[i].leftele - 1].ps.sldpth[k];
            totdpth += elem[river[i].leftele - 1].ps.sldpth[k];
        }
        for (k = 0; k < elem[river[i].rightele - 1].ps.nsoil; k++)
        {
            river[i].matl.bd += elem[river[i].rightele - 1].soil.bd[k] *
                elem[river[i].rightele - 1].ps.sldpth[k];
            totdpth += elem[river[i].rightele - 1].ps.sldpth[k];
        }

        river[i].matl.bd /= totdpth;
    }
}

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
            elem[i].np.no3 += MAX(elem[i].ns.no3[k], 0.0);
            elem[i].np.nh4 += MAX(elem[i].ns.nh4[k], 0.0);
        }

        NV_Ith(CV_Y, NO3(i)) = elem[i].np.no3;
        NV_Ith(CV_Y, NH4(i)) = elem[i].np.nh4;

        elem[i].ns0 = elem[i].ns;
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].ns.streamno3 = MAX(river[i].restart_input.streamno3, 0.0);
        river[i].ns.streamnh4 = MAX(river[i].restart_input.streamnh4, 0.0);
        river[i].ns.bedno3 = MAX(river[i].restart_input.bedno3, 0.0);
        river[i].ns.bednh4 = MAX(river[i].restart_input.bednh4, 0.0);

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
//    restart_file = fopen(restart_fn, "wb");
//    CheckFile(restart_file, restart_fn);
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
