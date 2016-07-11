#include "pihm.h"

void InitBGC (elem_struct *elem, int numele, river_struct *riv, int numriv, const epctbl_struct *epctbl, const ctrl_struct *ctrl)
{
    int             i;
    int             epc_ind;

    /* Detect if model is running in ensemble mode */
    if (verbose_mode)
    {
        printf ("BGC: Initializing BGC structures\n");
    }

    for (i = 0; i < numele; i++)
    {
        epc_ind = elem[i].attrib.lc_type - 1;

        if (epc_ind != ENF - 1 &&
            epc_ind != EBF - 1 &&
            epc_ind != DNF - 1 &&
            epc_ind != DBF - 1 &&
            epc_ind != GRASS - 1 &&
            epc_ind != CLOSE_SHRUB - 1 &&
            epc_ind != OPEN_SHRUB - 1)
        {
            printf ("Land cover type %d not been defined in Flux-PIHM-BGC\n",
                elem[i].attrib.lc_type);
            PihmExit (1);
        }

        elem[i].epc.woody = epctbl->woody[epc_ind];
        elem[i].epc.evergreen = epctbl->evergreen[epc_ind];
        elem[i].epc.c3_flag = epctbl->c3_flag[epc_ind];
        elem[i].epc.phenology_flag = epctbl->phenology_flag[epc_ind];
        elem[i].epc.onday = epctbl->onday[epc_ind];
        elem[i].epc.offday = epctbl->offday[epc_ind];
        elem[i].epc.transfer_days = epctbl->transfer_days[epc_ind];
        elem[i].epc.litfall_days = epctbl->litfall_days[epc_ind];
        elem[i].epc.leaf_turnover = epctbl->leaf_turnover[epc_ind];
        elem[i].epc.froot_turnover = epctbl->froot_turnover[epc_ind];
        elem[i].epc.livewood_turnover = epctbl->livewood_turnover[epc_ind];
        elem[i].epc.daily_mortality_turnover = epctbl->daily_mortality_turnover[epc_ind];
        elem[i].epc.daily_fire_turnover = epctbl->daily_fire_turnover[epc_ind];
        elem[i].epc.alloc_frootc_leafc = epctbl->alloc_frootc_leafc[epc_ind];
        elem[i].epc.alloc_newstemc_newleafc = epctbl->alloc_newstemc_newleafc[epc_ind];
        elem[i].epc.alloc_newlivewoodc_newwoodc = epctbl->alloc_newlivewoodc_newwoodc[epc_ind];
        elem[i].epc.alloc_crootc_stemc = epctbl->alloc_crootc_stemc[epc_ind];
        elem[i].epc.alloc_prop_curgrowth = epctbl->alloc_prop_curgrowth[epc_ind];
        elem[i].epc.avg_proj_sla = epctbl->avg_proj_sla[epc_ind];
        elem[i].epc.sla_ratio = epctbl->sla_ratio[epc_ind];
        elem[i].epc.lai_ratio = epctbl->lai_ratio[epc_ind];
        elem[i].epc.ext_coef = epctbl->ext_coef[epc_ind];
        elem[i].epc.flnr = epctbl->flnr[epc_ind];
        elem[i].epc.psi_open = epctbl->psi_open[epc_ind];
        elem[i].epc.psi_close = epctbl->psi_close[epc_ind];
        elem[i].epc.vpd_open = epctbl->vpd_open[epc_ind];
        elem[i].epc.vpd_close = epctbl->vpd_close[epc_ind];
        elem[i].epc.froot_cn = epctbl->froot_cn[epc_ind];
        elem[i].epc.leaf_cn = epctbl->leaf_cn[epc_ind];
        elem[i].epc.livewood_cn = epctbl->livewood_cn[epc_ind];
        elem[i].epc.deadwood_cn = epctbl->deadwood_cn[epc_ind];
        elem[i].epc.leaflitr_cn = epctbl->leaflitr_cn[epc_ind];
        elem[i].epc.leaflitr_flab = epctbl->leaflitr_flab[epc_ind];
        elem[i].epc.leaflitr_fucel = epctbl->leaflitr_fucel[epc_ind];
        elem[i].epc.leaflitr_fscel = epctbl->leaflitr_fscel[epc_ind];
        elem[i].epc.leaflitr_flig = epctbl->leaflitr_flig[epc_ind];
        elem[i].epc.frootlitr_flab = epctbl->frootlitr_flab[epc_ind];
        elem[i].epc.frootlitr_fucel = epctbl->frootlitr_fucel[epc_ind];
        elem[i].epc.frootlitr_fscel = epctbl->frootlitr_fscel[epc_ind];
        elem[i].epc.frootlitr_flig = epctbl->frootlitr_flig[epc_ind];
        elem[i].epc.deadwood_fucel = epctbl->deadwood_fucel[epc_ind];
        elem[i].epc.deadwood_fscel = epctbl->deadwood_fscel[epc_ind];
        elem[i].epc.deadwood_flig = epctbl->deadwood_flig[epc_ind];

        if (ctrl->spinup)
        {
            InitElemStor (&elem[i].stor, ctrl->spinupstart, ctrl->spinupend);
        }
    }

    for (i = 0; i < numriv; i++)
    {
        if (ctrl->spinup)
        {
            InitRiverStor (&riv[i].stor, ctrl->spinupstart, ctrl->spinupend);
        }
    }
}

void InitBGCVar (elem_struct *elem, int numele, river_struct *riv, int numriv, cinit_struct cinit, cstate_struct cs, nstate_struct ns, char *fn, int spinup)
{
    int             i;
    FILE           *init_file;

    /* Read initial conditions */
    if (!spinup)
    {
        init_file = fopen (fn, "rb");
        CheckFile (init_file, fn);

        for (i = 0; i < numele; i++)
        {
            fread (&elem[i].restart_input, sizeof (bgc_ic_struct), 1, init_file);

            restart_input (&elem[i].cs, &elem[i].ns, &elem[i].epv,
                &elem[i].restart_input);

            /* Calculate LAI for the coupling with Noah */
            elem[i].ps.proj_lai = elem[i].cs.leafc * elem[i].epc.avg_proj_sla;
        }

        for (i = 0; i < numriv; i++)
        {
            fread (&riv[i].sminn, sizeof (double), 1, init_file);
        }

        fclose (init_file);
    }
    else
    {
        for (i = 0; i < numele; i++)
        {
            elem[i].cinit = cinit;
            elem[i].cs = cs;
            elem[i].ns = ns;

            elem[i].ns.cwdn = cs.cwdc / elem[i].epc.deadwood_cn;
            elem[i].ns.litr2n = cs.litr2c / elem[i].epc.leaflitr_cn;
            elem[i].ns.litr3n = cs.litr3c / elem[i].epc.leaflitr_cn;
            elem[i].ns.litr4n = cs.litr4c / elem[i].epc.leaflitr_cn;
            elem[i].ns.soil1n = cs.soil1c / SOIL1_CN;
            elem[i].ns.soil2n = cs.soil2c / SOIL2_CN;
            elem[i].ns.soil3n = cs.soil3c / SOIL3_CN;
            elem[i].ns.soil4n = cs.soil4c / SOIL4_CN;

            if (elem[i].epc.evergreen)
            {
                elem[i].epv.dormant_flag = 0.0;
            }
            else
            {
                elem[i].epv.dormant_flag = 1.0;
            }

            //elem[i].epv.days_active = 0.;
            elem[i].epv.onset_flag = 0.0;
            elem[i].epv.onset_counter = 0.0;
            elem[i].epv.onset_gddflag = 0.0;
            elem[i].epv.onset_fdd = 0.0;
            elem[i].epv.onset_gdd = 0.0;
            elem[i].epv.onset_swi = 0.0;
            elem[i].epv.offset_flag = 0.0;
            elem[i].epv.offset_counter = 0.0;
            elem[i].epv.offset_fdd = 0.0;
            elem[i].epv.offset_swi = 0.0;
            elem[i].epv.annavg_t2m = elem[i].ps.tbot;

            firstday (&elem[i].epc, &elem[i].cinit, &elem[i].epv, &elem[i].cs, &elem[i].ns);
        }

        for (i = 0; i < numriv; i++)
        {
            riv[i].sminn = 0.0;
        }
    }

    for (i = 0; i < numele; i++)
    {
        zero_srcsnk (&elem[i].cs, &elem[i].ns, &elem[i].summary);
    }

    for (i = 0; i < numriv; i++)
    {
       riv[i].nleached_snk = 0.0;
       riv[i].sminn_leached = 0.0;
    }
}

//void Bgc2Noah (int t, pihm_struct pihm, lsm_struct noah, bgc_struct bgc)
//{
//    int             i;
//
//    for (i = 0; i < pihm->numele; i++)
//    {
//        if (!bgc->ctrl.spinup)
//        {
//            noah->grid[i].xlai = bgc->grid[i].epv.proj_lai;
//        }
//        else
//        {
//            if (pihm->elem[i].forc.lai_type > 0)
//            {
//                noah->grid[i].xlai = *pihm->elem[i].forc.lai;
//            }
//            else
//            {
//                noah->grid[i].xlai = MonthlyLAI (t, pihm->elem[i].lc.type);
//            }
//        }
//    }
//}

//void BgcCoupling (int t, int start_time, pihm_struct pihm, lsm_struct noah, bgc_struct bgc)
//{
//    static int      counter[MAXGRID];
//    static int      daylight_counter[MAXGRID];
//    double          dayl, prev_dayl;
//    spa_data        spa;
//    int             spa_result;
//    time_t          rawtime;
//    struct tm      *timestamp;
//    double          sfctmp;
//    double          solar;
//    int             i, k;
//    double          dummy[MAXGRID];
//    metvar_struct  *metv;
//    static int      first_balance;
//
//    if (t == start_time)
//    {
//        for (i = 0; i < pihm->numele; i++)
//        {
//            counter[i] = 0;
//            daylight_counter[i] = 0;
//
//            metv = &(bgc->grid[i].metv);
//            metv->tmax = -999.0;
//            metv->tmin = 999.0;
//            metv->tavg = 0.0;
//            metv->tsoil = 0.0;
//            metv->swc = 0.0;
//            metv->soilw = 0.0;
//            for (k = 0; k < 3; k++)
//            {
//                metv->latflux[k] = 0.0;
//            }
//            metv->tday = 0.0;
//            metv->q2d = 0.0;
//            metv->pa = 0.0;
//            metv->swavgfd = 0.0;
//            metv->par = 0.0;
//            metv->tnight = 0.0;
//        }
//        first_balance = 1;
//    }
//
//    for (i = 0; i < pihm->numele; i++)
//    {
//        metv = &(bgc->grid[i].metv);
//
//        sfctmp = noah->grid[i].sfctmp - 273.15;
//        metv->tmax = (metv->tmax > sfctmp) ? metv->tmax : sfctmp;
//        metv->tmin = (metv->tmin < sfctmp) ? metv->tmin : sfctmp;
//        metv->tavg += sfctmp;
//        solar = noah->grid[i].soldn;
//        metv->tsoil += noah->grid[i].stc[0] - 273.15;
//        metv->swc += noah->grid[i].soilw;
//        metv->soilw += noah->grid[i].soilm;
//        for (k = 0; k < 3; k++)
//        {
//            metv->latflux[k] += noah->grid[i].avgsubflux[k] * 1000.0 * 24.0 * 3600.0 / pihm->elem[i].topo.area;
//        }
//        metv->latflux[3] = 0.0;
//
//        if (solar > 1.0)
//        {
//            metv->tday += sfctmp;
//            metv->q2d += noah->grid[i].q2sat - noah->grid[i].q2;
//            metv->pa += noah->grid[i].sfcprs;
//            metv->swavgfd += solar;
//            metv->par += solar * RAD2PAR;
//            daylight_counter[i]++;
//        }
//        else
//        {
//            metv->tnight += sfctmp;
//        }
//
//        (counter[i])++;
//    }
//
//    if ((t - start_time) % 86400 == 0 && t > start_time)
//    {
//        rawtime = (int) (t - 86400);
//        timestamp = gmtime (&rawtime);
//        spa.year = timestamp->tm_year + 1900;
//        spa.month = timestamp->tm_mon + 1;
//        spa.day = timestamp->tm_mday;
//        spa.hour = timestamp->tm_hour;
//        spa.minute = timestamp->tm_min;
//        spa.second = timestamp->tm_sec;
//
//        spa.timezone = 0;
//        spa.delta_t = 67;
//        spa.delta_ut1 = 0;
//        spa.atmos_refract = 0.5667;
//
//        spa.longitude = bgc->grid[0].sitec.lon;
//        spa.latitude = bgc->grid[0].sitec.lat;
//        spa.elevation = 0.;
//        for (i = 0; i < pihm->numele; i++)
//        {
//            spa.elevation = spa.elevation + (double)pihm->elem[i].topo.zmax;
//        }
//        spa.elevation = spa.elevation / (double)pihm->numele;
//        /*
//         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
//         */
//        spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293., 5.26);
//        spa.temperature = noah->genprmt.tbot_data;
//
//        spa.function = SPA_ZA_RTS;
//        spa_result = spa_calculate (&spa);
//
//        /* daylength (s) */
//        dayl = (spa.sunset - spa.sunrise) * 3600.0;
//        dayl = (dayl < 0.0) ? (dayl + 24.0 * 3600.0) : dayl;
//
//        rawtime = rawtime - 24 * 3600;
//        timestamp = gmtime (&rawtime);
//        spa.year = timestamp->tm_year + 1900;
//        spa.month = timestamp->tm_mon + 1;
//        spa.day = timestamp->tm_mday;
//        spa.hour = timestamp->tm_hour;
//        spa.minute = timestamp->tm_min;
//        spa.second = timestamp->tm_sec;
//        spa_result = spa_calculate (&spa);
//        prev_dayl = (spa.sunset - spa.sunrise) * 3600.;
//        prev_dayl = (prev_dayl < 0.0) ? (prev_dayl + 12.0 * 3600.0) : prev_dayl;
//
//        for (i = 0; i < pihm->numele; i++)
//        {
//            metv = &(bgc->grid[i].metv);
//
//            metv->dayl = dayl;
//            metv->prev_dayl = prev_dayl;
//
//            metv->tavg /= (double) counter[i];
//            metv->tsoil /= (double) counter[i];
//            metv->swc /= (double) counter[i];
//            metv->soilw /= (double) counter[i];
//            for (k = 0; k < 3; k++)
//            {
//                metv->latflux[k] /= (double) counter[i];
//            }
//
//            metv->tday /= (double) daylight_counter[i];
//            metv->q2d /= (double) daylight_counter[i];
//            metv->pa /= (double) daylight_counter[i];
//            metv->swavgfd /= (double) daylight_counter[i];
//            metv->par /= (double) daylight_counter[i];
//
//            metv->tnight /= (double) (counter[i] - daylight_counter[i]);
//        }
//
//        DailyBgc (bgc, pihm->numele, pihm->numriv, t, dummy, first_balance);
//        first_balance = 0;
//
//        PrintData (bgc->prtctrl, bgc->ctrl.nprint, t, t - start_time, 86400,
//            pihm->ctrl.ascii);
//
//        
//        for (i = 0; i < pihm->numele; i++)
//        {
//            noah->grid[i].xlai = bgc->grid[i].epv.proj_lai;
//            noah->grid[i].cmcmax = noah->grid[i].cmcfactr * noah->grid[i].xlai;
//
//            metv = &(bgc->grid[i].metv);
//
//            counter[i] = 0;
//            daylight_counter[i] = 0;
//
//            metv->tmax = -999.0;
//            metv->tmin = 999.0;
//            metv->tavg = 0.0;
//            metv->tsoil = 0.0;
//            metv->swc = 0.0;
//            metv->soilw = 0.0;
//            for (k = 0; k < 3; k++)
//                metv->latflux[k] = 0.0;
//            metv->tday = 0.0;
//            metv->q2d = 0.0;
//            metv->pa = 0.0;
//            metv->swavgfd = 0.0;
//            metv->par = 0.0;
//            metv->tnight = 0.0;
//        }
//    }
//}

