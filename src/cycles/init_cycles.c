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
        int             kz, kcrop;
        /*
         * Initialize initial soil variables
         */
        soil_ind = elem[i].attrib.soil - 1;

        for (kz = 0; kz < MAXLYR; kz++)
        {
            if (kz < elem[i].ps.nlayers)
            {
                elem[i].soil.clay[kz] = soiltbl->clay_layer[soil_ind][kz];
                elem[i].soil.sand[kz] = soiltbl->sand_layer[soil_ind][kz];
                elem[i].soil.bd[kz]   = soiltbl->bd_layer[soil_ind][kz];
                elem[i].soil.fc[kz]   = soiltbl->fc[soil_ind][kz];
                elem[i].soil.pwp[kz]  = soiltbl->pwp[soil_ind][kz];
                elem[i].soil.b[kz]    = soiltbl->b[soil_ind][kz];
                elem[i].soil.air_entry_pot[kz] =
                                        soiltbl->air_entry_pot[soil_ind][kz];

                elem[i].soil.fc[kz]  *= cal->porosity;
                elem[i].soil.pwp[kz] *= cal->porosity;
            }
            else
            {
                elem[i].soil.clay[kz]          = BADVAL;
                elem[i].soil.sand[kz]          = BADVAL;
                elem[i].soil.bd[kz]            = BADVAL;
                elem[i].soil.fc[kz]            = BADVAL;
                elem[i].soil.pwp[kz]           = BADVAL;
                elem[i].soil.b[kz]             = BADVAL;
                elem[i].soil.air_entry_pot[kz] = BADVAL;
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
        int             kz;
        int             soil_ind;

        soil_ind = elem[i].attrib.soil - 1;

        elem[i].restart_input.water_residue_stan = 0.0;
        elem[i].restart_input.water_residue_flat = 0.0;
        elem[i].restart_input.c_residue_stan     = 0.0;
        elem[i].restart_input.c_residue_flat     = 0.0;
        elem[i].restart_input.c_manure_surface   = 0.0;
        elem[i].restart_input.n_residue_stan     = 0.0;
        elem[i].restart_input.n_residue_flat     = 0.0;
        elem[i].restart_input.n_manure_surface   = 0.0;

        for (kz = 0; kz < MAXLYR; kz++)
        {
            elem[i].restart_input.c_residue_abgd[kz]  = BADVAL;
            elem[i].restart_input.c_residue_root[kz]  = BADVAL;
            elem[i].restart_input.c_residue_rhizo[kz] = BADVAL;
            elem[i].restart_input.c_manure[kz]        = BADVAL;
            elem[i].restart_input.n_residue_abgd[kz]  = BADVAL;
            elem[i].restart_input.n_residue_root[kz]  = BADVAL;
            elem[i].restart_input.n_residue_rhizo[kz] = BADVAL;
            elem[i].restart_input.n_manure[kz]        = BADVAL;
            elem[i].restart_input.soc[kz]             = BADVAL;
            elem[i].restart_input.mbc[kz]             = BADVAL;
            elem[i].restart_input.son[kz]             = BADVAL;
            elem[i].restart_input.mbn[kz]             = BADVAL;
            elem[i].restart_input.no3[kz]             = BADVAL;
            elem[i].restart_input.nh4[kz]             = BADVAL;
        }

        for (kz = 0; kz < elem[i].ps.nlayers; kz++)
        {
            elem[i].restart_input.c_residue_abgd[kz]  = 0.0;
            elem[i].restart_input.c_residue_root[kz]  = 0.0;
            elem[i].restart_input.c_residue_rhizo[kz] = 0.0;
            elem[i].restart_input.c_manure[kz]        = 0.0;
            elem[i].restart_input.n_residue_abgd[kz]  = 0.0;
            elem[i].restart_input.n_residue_root[kz]  = 0.0;
            elem[i].restart_input.n_residue_rhizo[kz] = 0.0;
            elem[i].restart_input.n_manure[kz]        = 0.0;
            elem[i].restart_input.soc[kz]             = BADVAL;
            elem[i].restart_input.mbc[kz]             = BADVAL;
            elem[i].restart_input.son[kz]             = BADVAL;
            elem[i].restart_input.mbn[kz]             = BADVAL;
            elem[i].restart_input.no3[kz]             = BADVAL;
            elem[i].restart_input.nh4[kz]             = BADVAL;

            elem[i].restart_input.soc[kz] = soiltbl->om_layer[soil_ind][kz] /
                100.0 * 0.58 * elem[i].ps.soil_depth[kz] * elem[i].soil.bd[kz] *
                1.0E4;
            /* Initializes as 3% of SOC_Mass but "added" C */
            elem[i].restart_input.mbc[kz] = 0.03 *
                elem[i].restart_input.soc[kz];
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.son[kz] = elem[i].restart_input.soc[kz] * 0.1;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.mbn[kz] = elem[i].restart_input.mbc[kz] * 0.1;
            elem[i].restart_input.no3[kz] = soiltbl->no3[soil_ind][kz] *
                1.0E-3 * elem[i].ps.soil_depth[kz] / ctrl->soil_depth[kz];
            elem[i].restart_input.nh4[kz] = soiltbl->nh4[soil_ind][kz] *
                1.0E-3 * elem[i].ps.soil_depth[kz] / ctrl->soil_depth[kz];
        }
    }
}

void InitAgVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        int             kz;

        elem[i].ws.stan_residue   = elem[i].restart_input.water_residue_stan;
        elem[i].ws.flat_residue   = elem[i].restart_input.water_residue_flat;
        elem[i].cs.stan_residue   = elem[i].restart_input.c_residue_stan;
        elem[i].cs.flat_residue   = elem[i].restart_input.c_residue_flat;
        elem[i].cs.manure_surface = elem[i].restart_input.c_manure_surface;

        for (kz = 0; kz < MAXLYR; kz++)
        {
            elem[i].cs.abgd_residue[kz] =
                elem[i].restart_input.c_residue_abgd[kz];
            elem[i].cs.root_residue[kz] =
                elem[i].restart_input.c_residue_root[kz];
            elem[i].cs.rhizo_residue[kz] =
                elem[i].restart_input.c_residue_rhizo[kz];
            elem[i].cs.manure[kz] = elem[i].restart_input.c_manure[kz];
            elem[i].ns.residue_abgd[kz] =
                elem[i].restart_input.n_residue_abgd[kz];
            elem[i].ns.residue_root[kz] =
                elem[i].restart_input.n_residue_root[kz];
            elem[i].ns.residue_rhizo[kz] =
                elem[i].restart_input.n_residue_rhizo[kz];
            elem[i].ns.manure[kz] = elem[i].restart_input.n_manure[kz];
            elem[i].cs.soc[kz] = elem[i].restart_input.soc[kz];
            elem[i].cs.mbc[kz] = elem[i].restart_input.mbc[kz];
            elem[i].ns.son[kz] = elem[i].restart_input.son[kz];
            elem[i].ns.mbn[kz] = elem[i].restart_input.mbn[kz];
            elem[i].ns.no3[kz] = elem[i].restart_input.no3[kz];
            elem[i].ns.nh4[kz] = elem[i].restart_input.nh4[kz];
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