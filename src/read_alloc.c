#include "pihm.h"

void ReadAlloc(pihm_struct pihm)
{
    char            proj[MAXSTRING];
    char           *token;

    PIHMprintf(VL_VERBOSE, "\nRead input files:\n");

    strcpy(proj, project);
    if (strstr(proj, ".") != 0)
    {
        token = strtok(proj, ".");
        strcpy(proj, token);
    }
    else
    {
        strcpy(proj, project);
    }

    /* Set file names of the input files */
    sprintf(pihm->filename.riv,      "input/%s/%s.riv",      proj, proj);
    sprintf(pihm->filename.mesh,     "input/%s/%s.mesh",     proj, proj);
    sprintf(pihm->filename.att,      "input/%s/%s.att",      proj, proj);
    sprintf(pihm->filename.soil,     "input/%s/%s.soil",     proj, proj);
    sprintf(pihm->filename.lc,       "input/vegprmt.tbl");
    sprintf(pihm->filename.meteo,    "input/%s/%s.meteo",    proj, proj);
    sprintf(pihm->filename.lai,      "input/%s/%s.lai",      proj, proj);
    sprintf(pihm->filename.bc,       "input/%s/%s.bc",       proj, proj);
    sprintf(pihm->filename.para,     "input/%s/%s.para",     proj, proj);
    sprintf(pihm->filename.calib,    "input/%s/%s.calib",    proj, project);
    sprintf(pihm->filename.ic,       "input/%s/%s.ic",       proj, proj);
    sprintf(pihm->filename.tecplot,  "input/%s/%s.tecplot",  proj, proj);
#if defined(_FBR_)
    sprintf(pihm->filename.geol,     "input/%s/%s.geol",     proj, proj);
    sprintf(pihm->filename.bedrock,  "input/%s/%s.bedrock",  proj, proj);
#endif
#if defined(_NOAH_)
    sprintf(pihm->filename.lsm,      "input/%s/%s.lsm",      proj, proj);
    sprintf(pihm->filename.rad,      "input/%s/%s.rad",      proj, proj);
#endif
#if defined(_CYCLES_)
    sprintf(pihm->filename.cycles,   "input/%s/%s.cycles",   proj, proj);
    sprintf(pihm->filename.soilinit, "input/%s/%s.soilinit", proj, proj);
    sprintf(pihm->filename.crop,     "input/%s/%s.crop",     proj, proj);
    sprintf(pihm->filename.cyclesic, "input/%s/%s.cyclesic", proj, proj);
#endif
#if defined(_BGC_)
    sprintf(pihm->filename.bgc,      "input/%s/%s.bgc",      proj, proj);
    sprintf(pihm->filename.bgcic,    "input/%s/%s.bgcic",    proj, proj);
#endif

    /* Read river input file */
    ReadRiver(pihm->filename.riv, &pihm->rivtbl, &pihm->shptbl, &pihm->matltbl,
        &pihm->forc);

    /* Read mesh structure input file */
    ReadMesh(pihm->filename.mesh, &pihm->meshtbl);

    /* Read attribute table input file */
    ReadAtt(pihm->filename.att, &pihm->atttbl);

    /* Read soil input file */
    ReadSoil(pihm->filename.soil, &pihm->soiltbl);

    /* Read land cover input file */
    ReadLc(pihm->filename.lc, &pihm->lctbl);

    /* Read meteorological forcing input file */
    ReadForc(pihm->filename.meteo, &pihm->forc);

    /* Read LAI input file */
    ReadLai(pihm->filename.lai, &pihm->forc, &pihm->atttbl);

    /* Read source and sink input file */
    pihm->forc.nsource = 0;
#if NOT_YET_IMPLEMENTED
    ReadSS ();
#endif

    /* Read model control file */
    ReadPara(pihm->filename.para, &pihm->ctrl);

    /* Read calibration input file */
    ReadCalib(pihm->filename.calib, &pihm->cal);

    if (tecplot)
    {
        ReadTecplot(pihm->filename.tecplot, &pihm->ctrl);
    }

#if defined(_FBR_)
    /* Read geology input file */
    ReadGeol (pihm->filename.geol, &pihm->geoltbl);

    /* Read bedrock control file */
    ReadBedrock(pihm->filename.bedrock, &pihm->atttbl, &pihm->meshtbl,
        &pihm->ctrl);
#endif

    /* Read boundary condition input file
     * Boundary conditions might be needed by fbr thus should be read in after
     * reading bedrock input */
    ReadBc(pihm->filename.bc, &pihm->forc, &pihm->atttbl);

#if defined(_NOAH_)
    /* Read LSM input file */
    ReadLsm(pihm->filename.lsm, &pihm->siteinfo, &pihm->ctrl, &pihm->noahtbl);

    if (pihm->ctrl.rad_mode == TOPO_SOL)
    {
        /* Read radiation input file */
        ReadRad(pihm->filename.rad, &pihm->forc);
    }
#endif

#if defined(_CYCLES_)
    /* Read Cycles simulation control file */
    ReadCyclesCtrl(pihm->filename.cycles, &pihm->agtbl, &pihm->ctrl);

    /* Read soil initialization file */
    ReadSoilInit(pihm->filename.soilinit, &pihm->soiltbl);

    /* Read crop description file */
    ReadCrop(pihm->filename.crop, &pihm->croptbl);

    /* Read operation file */
    ReadOperation(&pihm->agtbl, &pihm->mgmttbl, &pihm->croptbl);
#endif

#if defined(_BGC_)
    ReadBgc(pihm->filename.bgc, &pihm->ctrl, &pihm->co2, &pihm->ndepctrl,
        &pihm->cninit, pihm->filename.co2, pihm->filename.ndep);

    /* Read Biome-BGC epc files */
    ReadEpc(&pihm->epctbl);

    /* Read CO2 and Ndep files */
    pihm->forc.co2 = (tsdata_struct *)malloc(sizeof(tsdata_struct));
    pihm->forc.ndep = (tsdata_struct *)malloc(sizeof(tsdata_struct));

    if (pihm->co2.varco2)
    {
        ReadAnnFile(&pihm->forc.co2[0], pihm->filename.co2);
    }

    if (pihm->ndepctrl.varndep)
    {
        ReadAnnFile(&pihm->forc.ndep[0], pihm->filename.ndep);
    }
#endif
}

void FreeData(pihm_struct pihm)
{
    int             i, j;

    /* Free river input structure */
    free(pihm->rivtbl.fromnode);
    free(pihm->rivtbl.tonode);
    free(pihm->rivtbl.down);
    free(pihm->rivtbl.leftele);
    free(pihm->rivtbl.rightele);
    free(pihm->rivtbl.shp);
    free(pihm->rivtbl.matl);
    free(pihm->rivtbl.bc);
    free(pihm->rivtbl.rsvr);

    free(pihm->shptbl.depth);
    free(pihm->shptbl.intrpl_ord);
    free(pihm->shptbl.coeff);

    free(pihm->matltbl.rough);
    free(pihm->matltbl.cwr);
    free(pihm->matltbl.ksath);
    free(pihm->matltbl.ksatv);
    free(pihm->matltbl.bedthick);

    /* Free mesh input structure */
    for (i = 0; i < nelem; i++)
    {
        free(pihm->meshtbl.node[i]);
        free(pihm->meshtbl.nabr[i]);
    }
    free(pihm->meshtbl.node);
    free(pihm->meshtbl.nabr);
    free(pihm->meshtbl.x);
    free(pihm->meshtbl.y);
    free(pihm->meshtbl.zmin);
    free(pihm->meshtbl.zmax);
#if defined(_FBR_)
    free(pihm->meshtbl.zbed);
#endif

    /* Free attribute input structure */
    for (i = 0; i < nelem; i++)
    {
        free(pihm->atttbl.bc[i]);
#if defined(_FBR_)
        free(pihm->atttbl.fbr_bc[i]);
#endif
    }
    free(pihm->atttbl.bc);
#if defined(_FBR_)
    free(pihm->atttbl.fbr_bc);
#endif
    free(pihm->atttbl.soil);
    free(pihm->atttbl.geol);
    free(pihm->atttbl.lc);
    free(pihm->atttbl.meteo);
    free(pihm->atttbl.lai);
    free(pihm->atttbl.source);

    /* Free soil input structure */
    free(pihm->soiltbl.silt);
    free(pihm->soiltbl.clay);
    free(pihm->soiltbl.om);
    free(pihm->soiltbl.bd);
    free(pihm->soiltbl.kinfv);
    free(pihm->soiltbl.ksatv);
    free(pihm->soiltbl.ksath);
    free(pihm->soiltbl.smcmax);
    free(pihm->soiltbl.smcmin);
    free(pihm->soiltbl.qtz);
    free(pihm->soiltbl.alpha);
    free(pihm->soiltbl.beta);
    free(pihm->soiltbl.areafh);
    free(pihm->soiltbl.areafv);
    free(pihm->soiltbl.dmac);
    free(pihm->soiltbl.smcref);
    free(pihm->soiltbl.smcwlt);

#if defined(_FBR_)
    /* Free geol input structure */
    free (pihm->geoltbl.ksatv);
    free (pihm->geoltbl.ksath);
    free (pihm->geoltbl.smcmax);
    free (pihm->geoltbl.smcmin);
    free (pihm->geoltbl.alpha);
    free (pihm->geoltbl.beta);
#endif

    /* Free landcover input structure */
    free(pihm->lctbl.laimax);
    free(pihm->lctbl.laimin);
    free(pihm->lctbl.vegfrac);
    free(pihm->lctbl.albedomin);
    free(pihm->lctbl.albedomax);
    free(pihm->lctbl.emissmin);
    free(pihm->lctbl.emissmax);
    free(pihm->lctbl.z0min);
    free(pihm->lctbl.z0max);
    free(pihm->lctbl.hs);
    free(pihm->lctbl.snup);
    free(pihm->lctbl.rgl);
    free(pihm->lctbl.rsmin);
    free(pihm->lctbl.rough);
    free(pihm->lctbl.rzd);

    /* Free forcing input structure */
    if (pihm->forc.nriverbc > 0)
    {
        for (i = 0; i < pihm->forc.nriverbc; i++)
        {
            for (j = 0; j < pihm->forc.riverbc[i].length; j++)
            {
                free(pihm->forc.riverbc[i].data[j]);
            }
            free(pihm->forc.riverbc[i].ftime);
            free(pihm->forc.riverbc[i].data);
        }
        free(pihm->forc.riverbc);
    }

    if (pihm->forc.nmeteo > 0)
    {
        for (i = 0; i < pihm->forc.nmeteo; i++)
        {
            for (j = 0; j < pihm->forc.meteo[i].length; j++)
            {
                free(pihm->forc.meteo[i].data[j]);
            }
            free(pihm->forc.meteo[i].ftime);
            free(pihm->forc.meteo[i].data);
            free(pihm->forc.meteo[i].value);
        }
        free(pihm->forc.meteo);
    }

    if (pihm->forc.nlai > 0)
    {
        for (i = 0; i < pihm->forc.nlai; i++)
        {
            for (j = 0; j < pihm->forc.lai[i].length; j++)
            {
                free(pihm->forc.lai[i].data[j]);
            }
            free(pihm->forc.lai[i].ftime);
            free(pihm->forc.lai[i].data);
            free(pihm->forc.lai[i].value);
        }
        free(pihm->forc.lai);
    }

    if (pihm->forc.nbc > 0)
    {
        for (i = 0; i < pihm->forc.nbc; i++)
        {
            for (j = 0; j < pihm->forc.bc[i].length; j++)
            {
                free(pihm->forc.bc[i].data[j]);
            }
            free(pihm->forc.bc[i].ftime);
            free(pihm->forc.bc[i].data);
            free(pihm->forc.bc[i].value);
        }
        free(pihm->forc.bc);
    }

#if defined(_NOAH_)
    if (pihm->forc.nrad > 0)
    {
        for (i = 0; i < pihm->forc.nrad; i++)
        {
            for (j = 0; j < pihm->forc.rad[i].length; j++)
            {
                free(pihm->forc.rad[i].data[j]);
            }
            free(pihm->forc.rad[i].ftime);
            free(pihm->forc.rad[i].data);
            free(pihm->forc.rad[i].value);
        }
        free(pihm->forc.rad);
    }
#endif

#if defined(_BGC_)
    free(pihm->epctbl.woody);
    free(pihm->epctbl.evergreen);
    free(pihm->epctbl.c3_flag);
    free(pihm->epctbl.phenology_flag);
    free(pihm->epctbl.onday);
    free(pihm->epctbl.offday);
    free(pihm->epctbl.transfer_days);
    free(pihm->epctbl.litfall_days);
    free(pihm->epctbl.leaf_turnover);
    free(pihm->epctbl.froot_turnover);
    free(pihm->epctbl.livewood_turnover);
    free(pihm->epctbl.daily_mortality_turnover);
    free(pihm->epctbl.daily_fire_turnover);
    free(pihm->epctbl.alloc_frootc_leafc);
    free(pihm->epctbl.alloc_newstemc_newleafc);
    free(pihm->epctbl.alloc_newlivewoodc_newwoodc);
    free(pihm->epctbl.alloc_crootc_stemc);
    free(pihm->epctbl.alloc_prop_curgrowth);
    free(pihm->epctbl.avg_proj_sla);
    free(pihm->epctbl.sla_ratio);
    free(pihm->epctbl.lai_ratio);
    free(pihm->epctbl.ext_coef);
    free(pihm->epctbl.flnr);
    free(pihm->epctbl.psi_open);
    free(pihm->epctbl.psi_close);
    free(pihm->epctbl.vpd_open);
    free(pihm->epctbl.vpd_close);
    free(pihm->epctbl.froot_cn);
    free(pihm->epctbl.leaf_cn);
    free(pihm->epctbl.livewood_cn);
    free(pihm->epctbl.deadwood_cn);
    free(pihm->epctbl.leaflitr_cn);
    free(pihm->epctbl.leaflitr_flab);
    free(pihm->epctbl.leaflitr_fucel);
    free(pihm->epctbl.leaflitr_fscel);
    free(pihm->epctbl.leaflitr_flig);
    free(pihm->epctbl.frootlitr_flab);
    free(pihm->epctbl.frootlitr_fucel);
    free(pihm->epctbl.frootlitr_fscel);
    free(pihm->epctbl.frootlitr_flig);
    free(pihm->epctbl.deadwood_fucel);
    free(pihm->epctbl.deadwood_fscel);
    free(pihm->epctbl.deadwood_flig);

    if (pihm->co2.varco2 > 0)
    {
        for (j = 0; j < pihm->forc.co2[0].length; j++)
        {
            free(pihm->forc.co2[0].data[j]);
        }
        free(pihm->forc.co2[0].ftime);
        free(pihm->forc.co2[0].data);
        free(pihm->forc.co2[0].value);
    }
    free(pihm->forc.co2);

    if (pihm->ndepctrl.varndep > 0)
    {
        for (j = 0; j < pihm->forc.ndep[0].length; j++)
        {
            free(pihm->forc.ndep[0].data[j]);
        }
        free(pihm->forc.ndep[0].ftime);
        free(pihm->forc.ndep[0].data);
        free(pihm->forc.ndep[0].value);
    }
    free(pihm->forc.ndep);
#endif

    free(pihm->ctrl.tout);

    /*
     * Close files
     */
    if (pihm->ctrl.waterbal)
    {
        fclose(pihm->print.watbal_file);
    }
    if (debug_mode)
    {
        fclose(pihm->print.cvodeperf_file);
    }
    for (i = 0; i < pihm->print.nprint; i++)
    {
        free(pihm->print.varctrl[i].var);
        free(pihm->print.varctrl[i].buffer);
        fclose(pihm->print.varctrl[i].datfile);
        if (pihm->ctrl.ascii)
        {
            fclose(pihm->print.varctrl[i].txtfile);
        }
    }
    if (tecplot)
    {
        for (i = 0; i < pihm->print.ntpprint; i++)
        {
            fclose(pihm->print.tp_varctrl[i].datfile);
        }
    }
    free(pihm->elem);
    free(pihm->river);
}
