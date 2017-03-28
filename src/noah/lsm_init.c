#include "pihm.h"

void InitLsm (elem_struct *elem, int numele, ctrl_struct ctrl,
    noahtbl_struct noahtbl, calib_struct cal)
{
    int             i;
    double          frzfact;

    for (i = 0; i < numele; i++)
    {

        /* Set-up soil layer depths */
        DefSldpth (elem[i].ps.sldpth, &elem[i].ps.nsoil, elem[i].ps.zsoil,
            elem[i].soil.depth, ctrl.sldpth, ctrl.nsoil);

        /* Set-up soil parameters */
        elem[i].ps.nmacd =
            FindLayer (elem[i].ps.sldpth, elem[i].ps.nsoil,
            elem[i].soil.dmac);

        elem[i].ps.nroot =
            FindLayer (elem[i].ps.sldpth, elem[i].ps.nsoil, elem[i].ps.rzd);

        RootDist (elem[i].ps.sldpth, elem[i].ps.nsoil, elem[i].ps.nroot,
            elem[i].ps.rtdis);
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
        frzfact =
            (elem[i].soil.smcmax / elem[i].soil.smcref) * (0.412 / 0.468);
        elem[i].ps.frzx = elem[i].ps.frzk * frzfact;
    }
}
