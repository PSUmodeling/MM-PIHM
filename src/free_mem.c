#include "pihm.h"

void FreeMem(pihm_struct pihm)
{
    int             i;

    FreeRivtbl(&pihm->rivtbl);

    FreeShptbl(&pihm->shptbl);

    FreeMatltbl(&pihm->matltbl);

    FreeMeshtbl(&pihm->meshtbl);

    FreeAtttbl(&pihm->atttbl);

    FreeSoiltbl(&pihm->soiltbl);

#if defined(_FBR_)
    FreeGeoltbl(&pihm->geoltbl);
#endif

    FreeLctbl(&pihm->lctbl);

    FreeForc(&pihm->forc);

#if defined(_BGC_)
    FreeEpctbl(&pihm->epctbl);
#endif

#if defined(_RT_)
    FreeRttbl(&pihm->rttbl);
#endif

    FreeCtrl(&pihm->ctrl);

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

void FreeRivtbl(rivtbl_struct *rivtbl)
{
    free(rivtbl->fromnode);
    free(rivtbl->tonode);
    free(rivtbl->down);
    free(rivtbl->leftele);
    free(rivtbl->rightele);
    free(rivtbl->shp);
    free(rivtbl->matl);
    free(rivtbl->bc);
    free(rivtbl->rsvr);
}

void FreeShptbl(shptbl_struct *shptbl)
{
    free(shptbl->depth);
    free(shptbl->intrpl_ord);
    free(shptbl->coeff);
}

void FreeMatltbl(matltbl_struct *matltbl)
{
    free(matltbl->rough);
    free(matltbl->cwr);
    free(matltbl->ksath);
    free(matltbl->ksatv);
    free(matltbl->bedthick);
}

void FreeMeshtbl(meshtbl_struct *meshtbl)
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        free(meshtbl->node[i]);
        free(meshtbl->nabr[i]);
    }
    free(meshtbl->node);
    free(meshtbl->nabr);
    free(meshtbl->x);
    free(meshtbl->y);
    free(meshtbl->zmin);
    free(meshtbl->zmax);
#if defined(_FBR_)
    free(meshtbl->zbed);
#endif
}

void FreeAtttbl(atttbl_struct *atttbl)
{
    int             i;

    /* Free attribute input structure */
    for (i = 0; i < nelem; i++)
    {
        free(atttbl->bc[i]);
#if defined(_FBR_)
        free(atttbl->fbr_bc[i]);
#endif
    }
    free(atttbl->bc);
#if defined(_FBR_)
    free(atttbl->fbr_bc);
#endif
    free(atttbl->soil);
    free(atttbl->geol);
    free(atttbl->lc);
    free(atttbl->meteo);
    free(atttbl->lai);
    free(atttbl->source);
}

void FreeSoiltbl(soiltbl_struct *soiltbl)
{
    free(soiltbl->silt);
    free(soiltbl->clay);
    free(soiltbl->om);
    free(soiltbl->bd);
    free(soiltbl->kinfv);
    free(soiltbl->ksatv);
    free(soiltbl->ksath);
    free(soiltbl->smcmax);
    free(soiltbl->smcmin);
    free(soiltbl->qtz);
    free(soiltbl->alpha);
    free(soiltbl->beta);
    free(soiltbl->areafh);
    free(soiltbl->areafv);
    free(soiltbl->dmac);
    free(soiltbl->smcref);
    free(soiltbl->smcwlt);
}

#if defined(_FBR_)
void FreeGeoltbl(geoltbl_struct *geoltbl)
{
    free(geoltbl->ksatv);
    free(geoltbl->ksath);
    free(geoltbl->smcmax);
    free(geoltbl->smcmin);
    free(geoltbl->alpha);
    free(geoltbl->beta);
}
#endif

void FreeLctbl(lctbl_struct *lctbl)
{
    /* Free landcover input structure */
    free(lctbl->laimax);
    free(lctbl->laimin);
    free(lctbl->vegfrac);
    free(lctbl->albedomin);
    free(lctbl->albedomax);
    free(lctbl->emissmin);
    free(lctbl->emissmax);
    free(lctbl->z0min);
    free(lctbl->z0max);
    free(lctbl->hs);
    free(lctbl->snup);
    free(lctbl->rgl);
    free(lctbl->rsmin);
    free(lctbl->rough);
    free(lctbl->rzd);
}

void FreeForc(forc_struct *forc)
{
    int             i, j;

    if (forc->nriverbc > 0)
    {
        for (i = 0; i < forc->nriverbc; i++)
        {
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                free(forc->riverbc[i].data[j]);
            }
            free(forc->riverbc[i].ftime);
            free(forc->riverbc[i].data);
        }
        free(forc->riverbc);
    }

    if (forc->nmeteo > 0)
    {
        for (i = 0; i < forc->nmeteo; i++)
        {
            for (j = 0; j < forc->meteo[i].length; j++)
            {
                free(forc->meteo[i].data[j]);
            }
            free(forc->meteo[i].ftime);
            free(forc->meteo[i].data);
            free(forc->meteo[i].value);
        }
        free(forc->meteo);
    }

    if (forc->nlai > 0)
    {
        for (i = 0; i < forc->nlai; i++)
        {
            for (j = 0; j < forc->lai[i].length; j++)
            {
                free(forc->lai[i].data[j]);
            }
            free(forc->lai[i].ftime);
            free(forc->lai[i].data);
            free(forc->lai[i].value);
        }
        free(forc->lai);
    }

    if (forc->nbc > 0)
    {
        for (i = 0; i < forc->nbc; i++)
        {
            for (j = 0; j < forc->bc[i].length; j++)
            {
                free(forc->bc[i].data[j]);
            }
            free(forc->bc[i].ftime);
            free(forc->bc[i].data);
            free(forc->bc[i].value);
        }
        free(forc->bc);
    }

#if defined(_NOAH_)
    if (forc->nrad > 0)
    {
        for (i = 0; i < forc->nrad; i++)
        {
            for (j = 0; j < forc->rad[i].length; j++)
            {
                free(forc->rad[i].data[j]);
            }
            free(forc->rad[i].ftime);
            free(forc->rad[i].data);
            free(forc->rad[i].value);
        }
        free(forc->rad);
    }
#endif

#if defined(_BGC_)
    if (forc->nco2 > 0)
    {
        for (j = 0; j < forc->co2[0].length; j++)
        {
            free(forc->co2[0].data[j]);
        }
        free(forc->co2[0].ftime);
        free(forc->co2[0].data);
        free(forc->co2[0].value);
    }
    free(forc->co2);

    if (forc->nndep > 0)
    {
        for (j = 0; j < forc->ndep[0].length; j++)
        {
            free(forc->ndep[0].data[j]);
        }
        free(forc->ndep[0].ftime);
        free(forc->ndep[0].data);
        free(forc->ndep[0].value);
    }
    free(forc->ndep);
#endif

#if defined(_RT_)
    if (forc->PrpFlg == 2)
    {
        for (i = 0; i < forc->nprcpc; i++)
        {
            for (j = 0; j < forc->prcpc[i].length; j++)
            {
                free(forc->prcpc[i].data[j]);
            }
            free(forc->prcpc[i].ftime);
            free(forc->prcpc[i].data);
            free(forc->prcpc[i].value);
        }
    }
#endif
}

#if defined(_BGC_)
void FreeEpctbl(epctbl_struct *epctbl)
{
    free(epctbl->woody);
    free(epctbl->evergreen);
    free(epctbl->c3_flag);
    free(epctbl->phenology_flag);
    free(epctbl->onday);
    free(epctbl->offday);
    free(epctbl->transfer_days);
    free(epctbl->litfall_days);
    free(epctbl->leaf_turnover);
    free(epctbl->froot_turnover);
    free(epctbl->livewood_turnover);
    free(epctbl->daily_mortality_turnover);
    free(epctbl->daily_fire_turnover);
    free(epctbl->alloc_frootc_leafc);
    free(epctbl->alloc_newstemc_newleafc);
    free(epctbl->alloc_newlivewoodc_newwoodc);
    free(epctbl->alloc_crootc_stemc);
    free(epctbl->alloc_prop_curgrowth);
    free(epctbl->avg_proj_sla);
    free(epctbl->sla_ratio);
    free(epctbl->lai_ratio);
    free(epctbl->ext_coef);
    free(epctbl->flnr);
    free(epctbl->psi_open);
    free(epctbl->psi_close);
    free(epctbl->vpd_open);
    free(epctbl->vpd_close);
    free(epctbl->froot_cn);
    free(epctbl->leaf_cn);
    free(epctbl->livewood_cn);
    free(epctbl->deadwood_cn);
    free(epctbl->leaflitr_cn);
    free(epctbl->leaflitr_flab);
    free(epctbl->leaflitr_fucel);
    free(epctbl->leaflitr_fscel);
    free(epctbl->leaflitr_flig);
    free(epctbl->frootlitr_flab);
    free(epctbl->frootlitr_fucel);
    free(epctbl->frootlitr_fscel);
    free(epctbl->frootlitr_flig);
    free(epctbl->deadwood_fucel);
    free(epctbl->deadwood_fscel);
    free(epctbl->deadwood_flig);
}
#endif

#if defined(_RT_)
void FreeRttbl(rttbl_struct *rttbl)
{
    if (rttbl->NumBTC > 0)
    {
        free(rttbl->BTC_loc);
    }

    if (rttbl->NumPUMP > 0)
    {
        free(rttbl->pumps);
    }
}
#endif

void FreeCtrl(ctrl_struct *ctrl)
{
    free(ctrl->tout);
}
