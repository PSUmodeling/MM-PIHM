#include "pihm.h"

void InitLsm(elem_struct *elem, const char ice_fn[], const ctrl_struct *ctrl,
    const noahtbl_struct *noahtbl, const calib_struct *cal)
{
    int             i;
    double          frzfact;
    int             read_ice_flag = 0;
    double         *iceh;
    FILE           *ice_file;
    const double    ICEH = 0.5;

    iceh = (double *)malloc(nelem * sizeof(double));

#if defined(_LUMPEDBGC_)
    for (i = 0; i < nelem + 1; i++)
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        if (elem[i].lc.glacier == 1)
        {
            read_ice_flag = 1;
            break;
        }
    }

    if (read_ice_flag == 1)
    {
        ice_file = fopen(ice_fn, "r");
        if (ice_file == NULL)
        {
            read_ice_flag = 0;
            pihm_printf(VL_NORMAL,
                "Optional input file *.ice is not available. "
                "Glacier ice depth will be initialized as 0.5 m.\n");
        }
        else
        {
            ReadGlacierIce(ice_fn, iceh);
        }
    }

#if defined(_LUMPEDBGC_)
    for (i = 0; i < nelem + 1; i++)
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        /* Set-up soil layer depths */
        DefSldpth(elem[i].ps.soil_depth, &elem[i].ps.nlayers, elem[i].ps.zsoil,
            elem[i].soil.depth, ctrl->soil_depth, ctrl->nlayers);

        /* Set-up glacier ice parameters */
        elem[i].ps.iceh = (elem[i].lc.glacier == 1) ?
            ((read_ice_flag == 1) ? iceh[i] : ICEH) : 0.0;

        /* Set-up soil parameters */
        elem[i].ps.nmacd =
            FindLayer(elem[i].ps.soil_depth, elem[i].ps.nlayers, elem[i].soil.dmac);

        elem[i].ps.nroot =
            FindLayer(elem[i].ps.soil_depth, elem[i].ps.nlayers, elem[i].ps.rzd);

        RootDist(elem[i].ps.soil_depth, elem[i].ps.nlayers, elem[i].ps.nroot,
            elem[i].ps.rtdis);

        /* Set-up universal parameters (not dependent on soil type or vegetation
         * type */
        elem[i].ps.sbeta = noahtbl->sbeta;
        elem[i].ps.salp = noahtbl->salp;
        elem[i].ps.frzk = noahtbl->frzk;
        elem[i].ps.fxexp = cal->fxexp * noahtbl->fxexp;
        elem[i].ps.czil = cal->czil * noahtbl->czil;
        elem[i].ps.lvcoef = noahtbl->lvcoef;
        elem[i].ps.zbot = noahtbl->zbot;
        elem[i].ps.tbot = noahtbl->tbot;

        /* To adjust frzk parameter to actual soil type */
        frzfact = (elem[i].soil.smcmax / elem[i].soil.smcref) * (0.412 / 0.468);
        elem[i].ps.frzx = elem[i].ps.frzk * frzfact;
    }

    free(iceh);
}
