#include "pihm.h"

void InitCycles(const calib_struct *cal, const agtbl_struct *agtbl,
    const mgmt_struct mgmttbl[], const crop_struct croptbl[],
    const soiltbl_struct *soiltbl, elem_struct elem[])
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
            if (k < elem[i].ps.nlayers)
            {
                elem[i].soil.clay[k] = soiltbl->clay_layer[soil_ind][k];
                elem[i].soil.sand[k] = soiltbl->sand_layer[soil_ind][k];
                elem[i].soil.bd[k]   = soiltbl->bd_layer[soil_ind][k];
                elem[i].soil.fc[k]   = soiltbl->fc[soil_ind][k];
                elem[i].soil.pwp[k]  = soiltbl->pwp[soil_ind][k];
                elem[i].soil.b[k]    = soiltbl->b[soil_ind][k];
                elem[i].soil.air_entry_pot[k] =
                                        soiltbl->air_entry_pot[soil_ind][k];

                elem[i].soil.fc[k]  *= cal->porosity;
                elem[i].soil.pwp[k] *= cal->porosity;
            }
            else
            {
                elem[i].soil.clay[k]          = BADVAL;
                elem[i].soil.sand[k]          = BADVAL;
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

void FirstDay(const soiltbl_struct *soiltbl, const ctrl_struct *ctrl,
    elem_struct elem[])
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        int             k;
        int             soil_ind;

        soil_ind = elem[i].attrib.soil_type - 1;

        elem[i].restart_input.water_residue_stan = 0.0;
        elem[i].restart_input.water_residue_flat = 0.0;
        elem[i].restart_input.c_residue_stan     = 0.0;
        elem[i].restart_input.c_residue_flat     = 0.0;
        elem[i].restart_input.c_manure_surface   = 0.0;
        elem[i].restart_input.n_residue_stan     = 0.0;
        elem[i].restart_input.n_residue_flat     = 0.0;
        elem[i].restart_input.n_manure_surface   = 0.0;

        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].restart_input.c_residue_abgd[k]  = BADVAL;
            elem[i].restart_input.c_residue_root[k]  = BADVAL;
            elem[i].restart_input.c_residue_rhizo[k] = BADVAL;
            elem[i].restart_input.c_manure[k]        = BADVAL;
            elem[i].restart_input.n_residue_abgd[k]  = BADVAL;
            elem[i].restart_input.n_residue_root[k]  = BADVAL;
            elem[i].restart_input.n_residue_rhizo[k] = BADVAL;
            elem[i].restart_input.n_manure[k]        = BADVAL;
            elem[i].restart_input.soc[k]             = BADVAL;
            elem[i].restart_input.mbc[k]             = BADVAL;
            elem[i].restart_input.son[k]             = BADVAL;
            elem[i].restart_input.mbn[k]             = BADVAL;
            elem[i].restart_input.no3[k]             = BADVAL;
            elem[i].restart_input.nh4[k]             = BADVAL;
        }

        for (k = 0; k < elem[i].ps.nlayers; k++)
        {
            elem[i].restart_input.c_residue_abgd[k]  = 0.0;
            elem[i].restart_input.c_residue_root[k]  = 0.0;
            elem[i].restart_input.c_residue_rhizo[k] = 0.0;
            elem[i].restart_input.c_manure[k]        = 0.0;
            elem[i].restart_input.n_residue_abgd[k]  = 0.0;
            elem[i].restart_input.n_residue_root[k]  = 0.0;
            elem[i].restart_input.n_residue_rhizo[k] = 0.0;
            elem[i].restart_input.n_manure[k]        = 0.0;
            elem[i].restart_input.soc[k]             = BADVAL;
            elem[i].restart_input.mbc[k]             = BADVAL;
            elem[i].restart_input.son[k]             = BADVAL;
            elem[i].restart_input.mbn[k]             = BADVAL;
            elem[i].restart_input.no3[k]             = BADVAL;
            elem[i].restart_input.nh4[k]             = BADVAL;

            elem[i].restart_input.soc[k] = soiltbl->om_layer[soil_ind][k] /
                100.0 * 0.58 * elem[i].ps.soil_depth[k] * elem[i].soil.bd[k] *
                1.0E4;
            /* Initializes as 3% of SOC_Mass but "added" C */
            elem[i].restart_input.mbc[k] = 0.03 * elem[i].restart_input.soc[k];
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.son[k] = elem[i].restart_input.soc[k] * 0.1;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.mbn[k] = elem[i].restart_input.mbc[k] * 0.1;
            elem[i].restart_input.no3[k] = soiltbl->no3[soil_ind][k] * 1.0E-3 *
                elem[i].ps.soil_depth[k] / ctrl->soil_depth[k];
            elem[i].restart_input.nh4[k] = soiltbl->nh4[soil_ind][k] * 1.0E-3 *
                elem[i].ps.soil_depth[k] / ctrl->soil_depth[k];
        }
    }
}

void InitAgVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        int             k;

        elem[i].ws.stan_residue   = elem[i].restart_input.water_residue_stan;
        elem[i].ws.flat_residue   = elem[i].restart_input.water_residue_flat;
        elem[i].cs.stan_residue   = elem[i].restart_input.c_residue_stan;
        elem[i].cs.flat_residue   = elem[i].restart_input.c_residue_flat;
        elem[i].cs.manure_surface = elem[i].restart_input.c_manure_surface;

        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].cs.abgd_residue[k] =
                elem[i].restart_input.c_residue_abgd[k];
            elem[i].cs.root_residue[k] =
                elem[i].restart_input.c_residue_root[k];
            elem[i].cs.rhizo_residue[k] =
                elem[i].restart_input.c_residue_rhizo[k];
            elem[i].cs.manure[k] = elem[i].restart_input.c_manure[k];
            elem[i].ns.residue_abgd[k] =
                elem[i].restart_input.n_residue_abgd[k];
            elem[i].ns.residue_root[k] =
                elem[i].restart_input.n_residue_root[k];
            elem[i].ns.residue_rhizo[k] =
                elem[i].restart_input.n_residue_rhizo[k];
            elem[i].ns.manure[k] = elem[i].restart_input.n_manure[k];
            elem[i].cs.soc[k] = elem[i].restart_input.soc[k];
            elem[i].cs.mbc[k] = elem[i].restart_input.mbc[k];
            elem[i].ns.son[k] = elem[i].restart_input.son[k];
            elem[i].ns.mbn[k] = elem[i].restart_input.mbn[k];
            elem[i].ns.no3[k] = elem[i].restart_input.no3[k];
            elem[i].ns.nh4[k] = elem[i].restart_input.nh4[k];
        }

        ZeroFluxes(&elem[i].wf, &elem[i].cf, &elem[i].nf);

        elem[i].ps.no3 = Profile(elem[i].ps.nlayers, elem[i].ns.no3);
        elem[i].ps.nh4 = Profile(elem[i].ps.nlayers, elem[i].ns.nh4);

        NV_Ith(CV_Y, SOLUTE_SOIL(i, NO3)) = elem[i].ps.no3;
        NV_Ith(CV_Y, SOLUTE_SOIL(i, NH4)) = elem[i].ps.nh4;

        elem[i].ns0 = elem[i].ns;
        elem[i].ps.no3_prev = elem[i].ps.no3;
        elem[i].ps.nh4_prev = elem[i].ps.nh4;
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].ns.no3 = 0.5 *
            (elem[river[i].left - 1].ps.no3 /
                elem[river[i].left - 1].soil.depth +
            elem[river[i].right - 1].ps.no3 /
                elem[river[i].right - 1].soil.depth) * river[i].ws.stage;

        river[i].ns.nh4 = 0.5 *
            (elem[river[i].left - 1].ps.nh4 /
                elem[river[i].left - 1].soil.depth +
            elem[river[i].right - 1].ps.nh4 /
                elem[river[i].right - 1].soil.depth) * river[i].ws.stage;

        NV_Ith(CV_Y, SOLUTE_RIVER(i, NO3)) = river[i].ns.no3;
        NV_Ith(CV_Y, SOLUTE_RIVER(i, NH4)) = river[i].ns.nh4;
    }
}
//
//void WriteCyclesIC(char *restart_fn, elem_struct *elem, river_struct *riv)
//{
//    int             i, j;
//    FILE           *restart_file;
//
//    restart_file = pihm_fopen(restart_fn, "wb");
//    pihm_printf(VL_VERBOSE, "Writing Cycles initial conditions.\n");
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