#include "pihm.h"

void InitLsm(const char ice_fn[], const ctrl_struct *ctrl, const noahtbl_struct *noahtbl, const calib_struct *calib,
    elem_struct elem[])
{
    int             i;
    double          frzfact;
    int             read_ice_flag = 0;
    double         *iceh;
    FILE           *fp;
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
        fp = fopen(ice_fn, "r");
        if (fp == NULL)
        {
            read_ice_flag = 0;
            pihm_printf(VL_NORMAL,
                "Optional input file *.ice is not available. Glacier ice depth will be initialized as 0.5 m.\n");
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
        // Set-up soil layer depths
        DefineSoilDepths(ctrl->nlayers, elem[i].soil.depth, ctrl->soil_depth, &elem[i].ps.nlayers,
            elem[i].ps.soil_depth, elem[i].ps.zsoil);

        // Set-up glacier ice parameters
        elem[i].ps.iceh = (elem[i].lc.glacier == 1) ? ((read_ice_flag == 1) ? iceh[i] : ICEH) : 0.0;

        // Set-up soil parameters
        elem[i].ps.nmacd = FindLayer(elem[i].ps.nlayers, elem[i].soil.dmac, elem[i].ps.soil_depth);

        elem[i].ps.nroot = FindLayer(elem[i].ps.nlayers, elem[i].ps.rzd, elem[i].ps.soil_depth);

        RootDist(elem[i].ps.nlayers, elem[i].ps.nroot, elem[i].ps.soil_depth, elem[i].ps.rtdis);

        // Set-up universal parameters (not dependent on soil type or vegetation type
        elem[i].ps.sbeta = noahtbl->sbeta;
        elem[i].ps.salp = noahtbl->salp;
        elem[i].ps.frzk = noahtbl->frzk;
        elem[i].ps.fxexp = calib->fxexp * noahtbl->fxexp;
        elem[i].ps.czil = calib->czil * noahtbl->czil;
        elem[i].ps.lvcoef = noahtbl->lvcoef;
        elem[i].ps.zbot = noahtbl->zbot;
        elem[i].ps.tbot = noahtbl->tbot;

        // To adjust frzk parameter to actual soil type
        frzfact = (elem[i].soil.smcmax / elem[i].soil.smcref) * 0.412 / 0.468;
        elem[i].ps.frzx = elem[i].ps.frzk * frzfact;
    }

    free(iceh);
}
