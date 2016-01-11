/*****************************************************************************
 * File		:   lsm_init.c 
 * Function	:   Noah initialization functions
 ****************************************************************************/
#include "pihm.h"

void InitLsm (elem_struct *elem, int numele, ctrl_struct ctrl,
    noahtbl_struct noahtbl, calib_struct cal)
{
    int             i;
    double          frzfact;

    for (i = 0; i < numele; i++)
    {

        /* Set-up soil layer depths */
        DefSldpth (elem[i].ps.sldpth, &elem[i].ps.nsoil, elem[i].soil.depth,
            ctrl.sldpth, ctrl.nsoil);

        /* Set-up soil parameters */
        elem[i].ps.nmacd =
            FindLayer (elem[i].ps.sldpth, elem[i].ps.nsoil,
            elem[i].soil.dmac);

        elem[i].lc.nroot = FindLayer (elem[i].ps.sldpth, elem[i].ps.nsoil, elem[i].lc.rzd);

        RootDist (elem[i].ps.sldpth, elem[i].ps.nsoil, elem[i].lc.nroot,
            elem[i].lc.rtdis);
        /* Set-up universal parameters (not dependent on soil type or
         * vegetation type */
        elem[i].ps.sbeta = noahtbl.sbeta;
        elem[i].ps.salp = noahtbl.salp;
        elem[i].ps.frzk = noahtbl.frzk;
        elem[i].ps.fxexp = cal.fxexp * noahtbl.fxexp;
        elem[i].ps.czil = cal.czil * noahtbl.czil;
        elem[i].ps.lvcoef = noahtbl.lvcoef;
        elem[i].ps.zbot = noahtbl.zbot;
        elem[i].ps.tbot = noahtbl.tbot;

        /* To adjust frzk parameter to actual soil type */
        frzfact = (elem[i].soil.smcmax / elem[i].soil.smcref) * (0.412 / 0.468);
        elem[i].ps.frzx = elem[i].ps.frzk * frzfact;
    }
}
//
////void LsmFreeData (pihm_struct pihm, lsm_struct noah)
////{
////    int             i, j;
////
////    if (noah->rad_mode == 1)
////    {
////        for (i = 0; i < noah->forcing.nts; i++)
////        {
////            for (j = 0; j < noah->forcing.ts[i].length; j++)
////            {
////                free (noah->forcing.ts[i].data[j]);
////            }
////            free (noah->forcing.ts[i].ftime);
////            free (noah->forcing.ts[i].data);
////        }
////        free (noah->forcing.ts);
////        for (i = 0; i < 2; i++)
////        {
////            free (noah->forcing.rad[i]);
////        }
////    }
////
////    free (noah->grid);
////    free (noah->ic.t1);
////    free (noah->ic.snowh);
////    for (i = 0; i < pihm->numele; i++)
////    {
////        free (noah->ic.stc[i]);
////        free (noah->ic.smc[i]);
////        free (noah->ic.sh2o[i]);
////    }
////    free (noah->ic.stc);
////    free (noah->ic.smc);
////    free (noah->ic.sh2o);
////
////    for (i = 0; i < noah->nprint; i++)
////    {
////        free (noah->prtctrl[i].vrbl);
////        free (noah->prtctrl[i].buffer);
////    }
////}
