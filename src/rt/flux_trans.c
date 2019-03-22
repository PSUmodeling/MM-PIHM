#include "pihm.h"

#define MIN(a,b) (((a)<(b))? (a):(b))
#define MAX(a,b) (((a)>(b))? (a):(b))

void fluxtrans(const pihm_struct pihm, Chem_Data CD)
{
    int             i, k = 0;
    int             BOUND_VOL = CD->NumVol;
    int             PRCP_VOL = CD->NumVol - 1;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          heqv;
        double          satn;
        int             k;

        UpdateVcele(MAX(pihm->elem[i].ws.gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].soil.smcmax,
            pihm->elem[i].soil.smcmin, pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        satn = UnsatSatRatio(pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_UNSAT(i)]);

#if defined(_FBR_)
        UpdateVcele(MAX(pihm->elem[i].ws.fbr_gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_FBR_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].geol.smcmax,
            pihm->elem[i].geol.smcmin, pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        satn = UnsatSatRatio(pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_FBR_UNSAT(i)]);
#endif

        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_GW(i)].chms.t_conc[k] = CD->Vcele[RT_GW(i)].chms.t_mole[k] /
                CD->Vcele[RT_GW(i)].vol / CD->Vcele[RT_GW(i)].porosity;

            CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] = CD->Vcele[RT_UNSAT(i)].chms.t_mole[k] /
                CD->Vcele[RT_UNSAT(i)].vol / CD->Vcele[RT_UNSAT(i)].porosity;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    /* Update river cells */
    for (i = 0; i < nriver; i++)
    {
        int             k;

        UpdateVcele(MAX(pihm->river[i].ws.gw, 1.0E-5) +
            MAX(pihm->river[i].ws.stage, 1.0E-5) /
            CD->Vcele[RT_RIVER(i)].porosity, 1.0, &CD->Vcele[RT_RIVER(i)]);

        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_RIVER(i)].chms.t_conc[k] = CD->Vcele[RT_RIVER(i)].chms.t_mole[k] /
                CD->Vcele[RT_RIVER(i)].vol / CD->Vcele[RT_RIVER(i)].porosity;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < 3; j++)
        {
            /* Flux for GW lateral flow */
            CD->Flux[RT_LAT_GW(i, j)].flux = pihm->elem[i].wf.subsurf[j];

            /* Flux for UNSAT lateral flow */
            CD->Flux[RT_LAT_UNSAT(i, j)].s_area = 0.5 *
                pihm->elem[i].topo.edge[j] *
                (CD->Vcele[CD->Flux[RT_LAT_UNSAT(i, j)].nodeup - 1].height_t +
                CD->Vcele[CD->Flux[RT_LAT_UNSAT(i, j)].nodelo - 1].height_t);

#if defined(_FBR_)
            /* Flux for deep lateral flow */
            CD->Flux[RT_LAT_FBR_GW(i, j)].flux = pihm->elem[i].wf.fbrflow[j];

            /* Flux for bedrock unsat lateral flow */
            CD->Flux[RT_LAT_FBR_UNSAT(i, j)].s_area = 0.5 *
                pihm->elem[i].topo.edge[j] *
                (CD->Vcele[CD->Flux[RT_LAT_FBR_UNSAT(i, j)].nodeup - 1].height_t +
                CD->Vcele[CD->Flux[RT_LAT_FBR_UNSAT(i, j)].nodelo - 1].height_t);

#endif
        }

        /* Flux for UNSAT - GW vertical flow */
        CD->Flux[RT_RECHG_UNSAT(i)].flux = pihm->elem[i].wf.rechg *
            CD->Vcele[RT_UNSAT(i)].area;

        CD->Flux[RT_RECHG_GW(i)].flux = -pihm->elem[i].wf.rechg *
            CD->Vcele[RT_GW(i)].area;

        CD->Flux[RT_INFIL(i)].flux =
            -((pihm->elem[i].wf.infil > 0.0) ? pihm->elem[i].wf.infil : 0.0) *
            pihm->elem[i].topo.area;

#if defined(_FBR_)
        CD->Flux[RT_FBR_RECHG_UNSAT(i)].flux = pihm->elem[i].wf.fbr_rechg * CD->Vcele[RT_FBR_RECHG_UNSAT(i)].area;
        CD->Flux[RT_FBR_RECHG_GW(i)].flux = -pihm->elem[i].wf.fbr_rechg * CD->Vcele[RT_FBR_RECHG_GW(i)].area;
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[RT_LEFT_SURF2RIVER(i)].flux = pihm->river[i].wf.rivflow[2];
        CD->Flux[RT_RIGHT_SURF2RIVER(i)].flux = pihm->river[i].wf.rivflow[3];
        CD->Flux[RT_LEFT_AQIF2RIVER(i)].flux = pihm->river[i].wf.rivflow[7] +
            pihm->river[i].wf.rivflow[4];
        CD->Flux[RT_RIGHT_AQIF2RIVER(i)].flux = pihm->river[i].wf.rivflow[8] +
            pihm->river[i].wf.rivflow[5];
        CD->Flux[RT_DOWN_RIVER2RIVER(i)].flux = pihm->river[i].wf.rivflow[9] +
            pihm->river[i].wf.rivflow[1];
        CD->Flux[RT_UP_RIVER2RIVER(i)].flux = pihm->river[i].wf.rivflow[10] +
            pihm->river[i].wf.rivflow[0];

        if (CD->Flux[RT_UP_RIVER2RIVER(i)].node_trib > 0)
        {
            CD->Flux[RT_UP_RIVER2RIVER(i)].flux_trib =
                -(pihm->river[pihm->river[i].up[1] - 1].wf.rivflow[9] +
                pihm->river[pihm->river[i].up[1] - 1].wf.rivflow[1]);
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        /* For gw cells, contact area is needed for dispersion; */
        for (j = 0; j < 3; j++)
        {
            double              h1, h2;

            if (CD->Flux[RT_LAT_GW(i, j)].BC == NO_FLOW)
            {
                continue;
            }

            if (pihm->elem[i].nabr[j] > 0)
            {
                h1 = 0.5 *
                    (CD->Vcele[RT_GW(i)].height_o +
                    CD->Vcele[RT_GW(i)].height_t);
                h2 = 0.5 *
                    (CD->Vcele[RT_GW(pihm->elem[i].nabr[j] - 1)].height_o +
                    CD->Vcele[RT_GW(pihm->elem[i].nabr[j] - 1)].height_t);

                CD->Flux[RT_LAT_GW(i, j)].s_area =
                    (CD->Flux[RT_LAT_GW(i, j)].flux > 0.0) ?
                    pihm->elem[i].topo.edge[j] * h1 :
                    pihm->elem[i].topo.edge[j] * h2;
            }
            else if (pihm->elem[i].nabr[j] < 0)
            {
                h1 = 0.5 *
                    (CD->Vcele[RT_GW(i)].height_o +
                    CD->Vcele[RT_GW(i)].height_t);
                h2 = 0.5 *
                    (CD->Vcele[RT_RIVER(-pihm->elem[i].nabr[j] - 1)].height_o +
                    CD->Vcele[RT_RIVER(-pihm->elem[i].nabr[j] - 1)].height_t);

                CD->Flux[RT_LAT_GW(i, j)].s_area =
                    (CD->Flux[RT_LAT_GW(i, j)].flux > 0.0) ?
                    pihm->elem[i].topo.edge[j] * h1 :
                    pihm->elem[i].topo.edge[j] * h2;
            }

            /* Calculate velocity according to flux and area */
            CD->Flux[RT_LAT_GW(i, j)].velocity =
                (CD->Flux[RT_LAT_GW(i, j)].s_area > 1.0E-4) ?
                CD->Flux[RT_LAT_GW(i, j)].flux /
                CD->Flux[RT_LAT_GW(i, j)].s_area : 1.0E-10 / 86400;
        }

        CD->Flux[RT_RECHG_UNSAT(i)].s_area = pihm->elem[i].topo.area;
        CD->Flux[RT_RECHG_UNSAT(i)].velocity =
            CD->Flux[RT_RECHG_UNSAT(i)].flux / pihm->elem[i].topo.area;

        CD->Flux[RT_RECHG_GW(i)].s_area = pihm->elem[i].topo.area;
        CD->Flux[RT_RECHG_GW(i)].velocity =
            CD->Flux[RT_RECHG_GW(i)].flux / pihm->elem[i].topo.area;
    }

    /* Correct river flux area and velocity */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].leftele - 1].nabr[j] == i + 1)
            {
                CD->Flux[RT_LEFT_AQIF2RIVER(i)].s_area =
                CD->Flux[RT_LAT_GW(pihm->river[i].leftele - 1, j)].s_area;
                CD->Flux[RT_LEFT_AQIF2RIVER(i)].velocity =
                -CD->Flux[RT_LAT_GW(pihm->river[i].leftele - 1, j)].velocity;
                break;
            }
        }

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].rightele - 1].nabr[j] == i + 1)
            {
                CD->Flux[RT_RIGHT_AQIF2RIVER(i)].s_area =
                CD->Flux[RT_LAT_GW(pihm->river[i].rightele - 1, j)].s_area;
                CD->Flux[RT_RIGHT_AQIF2RIVER(i)].velocity =
                -CD->Flux[RT_LAT_GW(pihm->river[i].rightele - 1, j)].velocity;
                break;
            }
        }
    }

    /* Update virtual volume */

    if (pihm->ctrl.PrpFlg)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[PRCP_VOL - 1].chms.t_conc[k] =
                pihm->rttbl.prcp_conc[k] * pihm->rttbl.Condensation;
        }
    }
    else
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[PRCP_VOL - 1].chms.t_conc[k] = 0.0;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < pihm->rttbl.NumStc; k++)
    {
        CD->Vcele[BOUND_VOL - 1].chms.t_conc[k] =
            pihm->rttbl.prcp_conc[k] * pihm->rttbl.Condensation;
        CD->Vcele[BOUND_VOL - 1].chms.p_conc[k] =
            pihm->rttbl.prcp_conc[k] * pihm->rttbl.Condensation;
    }

    /*
     * Transport
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < NumSpc; j++)
        {
            if (pihm->chemtbl[j].mtype == MIXED_MA)
            {
                for (k = 0; k < pihm->rttbl.NumSsc; k++)
                {
                    if ((pihm->rttbl.Totalconc[j][k + pihm->rttbl.NumStc] != 0) &&
                        (pihm->chemtbl[k + pihm->rttbl.NumStc].itype != AQUEOUS))
                    {
                        CD->Vcele[RT_GW(i)].chms.t_conc[j] -=
                            pihm->rttbl.Totalconc[j][k + pihm->rttbl.NumStc] *
                            CD->Vcele[RT_GW(i)].chms.s_conc[k];
                        CD->Vcele[RT_UNSAT(i)].chms.t_conc[j] -=
                            pihm->rttbl.Totalconc[j][k + pihm->rttbl.NumStc] *
                            CD->Vcele[RT_UNSAT(i)].chms.s_conc[k];
                    }
                }
            }
        }
    }

    OS3D(pihm->chemtbl, &pihm->rttbl, CD);
}

void SpeciationReaction(int t, int stepsize, const pihm_struct pihm,
    Chem_Data CD)
{
    int             i, k;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          heqv;
        double          satn;
        int             k;

        UpdateVcele(MAX(pihm->elem[i].ws.gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].soil.smcmax,
            pihm->elem[i].soil.smcmin, pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        satn = UnsatSatRatio(pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_UNSAT(i)]);

#if defined(_FBR_)
        UpdateVcele(MAX(pihm->elem[i].ws.fbr_gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_FBR_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].geol.smcmax,
            pihm->elem[i].geol.smcmin, pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        satn = UnsatSatRatio(pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_FBR_UNSAT(i)]);
#endif

        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_GW(i)].chms.t_conc[k] = CD->Vcele[RT_GW(i)].chms.t_mole[k] /
                CD->Vcele[RT_GW(i)].vol / CD->Vcele[RT_GW(i)].porosity;
            CD->Vcele[RT_GW(i)].chms.t_conc[k] =
                (CD->Vcele[RT_GW(i)].chms.t_conc[k] > 1.0E-20) ?
                CD->Vcele[RT_GW(i)].chms.t_conc[k] : 1.0E-20;

            CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] = CD->Vcele[RT_UNSAT(i)].chms.t_mole[k] /
                CD->Vcele[RT_UNSAT(i)].vol / CD->Vcele[RT_UNSAT(i)].porosity;
            CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] =
                (CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] > 1.0E-20) ?
                CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] : 1.0E-20;
        }

        //for (j = 0; j < NumSpc; j++)
        //{
        //    if (pihm->chemtbl[j].mtype == MIXED_MA)
        //    {
        //        for (k = 0; k < pihm->rttbl.NumSsc; k++)
        //        {
        //            if ((CD->Totalconc[j][k + pihm->rttbl.NumStc] != 0) &&
        //                (pihm->chemtbl[k + pihm->rttbl.NumStc].itype != AQUEOUS))
        //            {
        //                CD->Vcele[RT_GW(i)].chms.t_conc[j] =
        //                    CD->Vcele[RT_GW(i)].chms.t_conc[j] + CD->Totalconc[j][k +
        //                    pihm->rttbl.NumStc] * CD->Vcele[RT_GW(i)].chms.s_conc[k];
        //                CD->Vcele[RT_UNSAT(i)].chms.t_conc[j] =
        //                    CD->Vcele[RT_UNSAT(i)].chms.t_conc[j] + CD->Totalconc[j][k +
        //                    pihm->rttbl.NumStc] * CD->Vcele[RT_UNSAT(i)].chms.s_conc[k];
        //            }
        //        }
        //    }
        //}
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    /* Update river cells */
    for (i = 0; i < nriver; i++)
    {
        int             k;

        UpdateVcele(MAX(pihm->river[i].ws.gw, 1.0E-5) +
            MAX(pihm->river[i].ws.stage, 1.0E-5) /
            CD->Vcele[RT_RIVER(i)].porosity, 1.0, &CD->Vcele[RT_RIVER(i)]);

        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_RIVER(i)].chms.t_conc[k] = CD->Vcele[RT_RIVER(i)].chms.t_mole[k] /
                CD->Vcele[RT_RIVER(i)].vol / CD->Vcele[RT_RIVER(i)].porosity;
            CD->Vcele[RT_RIVER(i)].chms.t_conc[k] =
                (CD->Vcele[RT_RIVER(i)].chms.t_conc[k] > 1.0E-20) ?
                CD->Vcele[RT_RIVER(i)].chms.t_conc[k] : 1.0E-20;
        }
    }

    if (t - pihm->ctrl.starttime >= pihm->ctrl.RT_delay)
    {
        /*
         * Reaction
         */
        if ((!pihm->rttbl.RecFlg) && (t > pihm->ctrl.starttime) &&
            (t - pihm->ctrl.starttime) % pihm->ctrl.AvgScl == 0)
        {
#ifdef _OPENMP
# pragma omp parallel for
#endif
            for (i = 0; i < nelem; i++)
            {
                double          t_conc0[MAXSPS];
                int             k;

                for (k = 0; k < NumSpc; k++)
                {
                    t_conc0[k] = CD->Vcele[RT_GW(i)].chms.t_conc[k];
                }
                React((double)pihm->ctrl.AvgScl, pihm->chemtbl, pihm->kintbl,
                    &pihm->rttbl, CD->Vcele[RT_GW(i)].sat,
                    &CD->Vcele[RT_GW(i)].chms);

                for (k = 0; k < NumSpc; k++)
                {
                    CD->Vcele[RT_GW(i)].chmf.react_flux[k] =
                        (CD->Vcele[RT_GW(i)].chms.t_conc[k] - t_conc0[k]) *
                        CD->Vcele[RT_GW(i)].vol /
                        (double)pihm->ctrl.AvgScl;
                }

                for (k = 0; k < NumSpc; k++)
                {
                    t_conc0[k] = CD->Vcele[RT_UNSAT(i)].chms.t_conc[k];
                }
                React((double)pihm->ctrl.AvgScl, pihm->chemtbl, pihm->kintbl,
                    &pihm->rttbl, CD->Vcele[RT_UNSAT(i)].sat,
                    &CD->Vcele[RT_UNSAT(i)].chms);

                for (k = 0; k < NumSpc; k++)
                {
                    CD->Vcele[RT_UNSAT(i)].chmf.react_flux[k] =
                        (CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] - t_conc0[k]) *
                        CD->Vcele[RT_UNSAT(i)].vol / (double)pihm->ctrl.AvgScl;
                }
            }
        }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             j;

            for (j = 0; j < pihm->rttbl.NumStc; j++)
            {
                if (pihm->chemtbl[j].itype == MINERAL)
                {
                    /* Averaging mineral concentration to ensure mass
                     * conservation !! */
                    CD->Vcele[RT_GW(i)].chms.t_conc[j] =
                        (CD->Vcele[RT_GW(i)].chms.t_conc[j] *
                        CD->Vcele[RT_GW(i)].height_t +
                        CD->Vcele[RT_UNSAT(i)].chms.t_conc[j] *
                        (pihm->elem[i].soil.depth -
                        CD->Vcele[RT_GW(i)].height_t)) /
                        pihm->elem[i].soil.depth;
                    CD->Vcele[RT_UNSAT(i)].chms.t_conc[j] =
                        CD->Vcele[RT_GW(i)].chms.t_conc[j];
                    CD->Vcele[RT_GW(i)].chms.p_conc[j] =
                        CD->Vcele[RT_GW(i)].chms.t_conc[j];
                    CD->Vcele[RT_UNSAT(i)].chms.p_conc[j] =
                        CD->Vcele[RT_GW(i)].chms.t_conc[j];
                }
            }
        }
    } /* RT step control ends */

    /* Reset fluxes for next averaging stage */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < CD->NumFac; k++)
    {
        CD->Flux[k].velocity = 0.0;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;
        CD->Flux[k].s_area = 0.0;
    }

    /* Every hour */
    if ((t - pihm->ctrl.starttime) % 3600 == 0)
    {
        if (!pihm->rttbl.RecFlg)
        {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
            for (i = 0; i < pihm->rttbl.NumStc; i++)
            {
                int             j;

                for (j = 0; j < nriver; j++)
                {
                    CD->Vcele[RT_RIVER(j)].chms.p_conc[i] =
                        (pihm->chemtbl[i].itype == MINERAL) ?
                        CD->Vcele[RT_RIVER(j)].chms.t_conc[i] :
                        fabs(CD->Vcele[RT_RIVER(j)].chms.t_conc[i] * 0.1);
                }
            }
        }

        if (!pihm->rttbl.RecFlg)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < nriver; i++)
            {
                double          t_conc0[MAXSPS];
                int             k;

                for (k = 0; k < NumSpc; k++)
                {
                    t_conc0[k] = CD->Vcele[RT_RIVER(i)].chms.t_conc[k];
                }

                Speciation(pihm->chemtbl, &pihm->rttbl, 0,
                    &CD->Vcele[RT_RIVER(i)].chms);

                for (k = 0; k < NumSpc; k++)
                {
                    CD->Vcele[RT_RIVER(i)].chmf.react_flux[k] =
                        (CD->Vcele[RT_RIVER(i)].chms.t_conc[k] - t_conc0[k]) *
                        CD->Vcele[RT_RIVER(i)].vol / 3600.0;
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < CD->NumVol; i++)
            {
                if (CD->Vcele[i].type != VIRTUAL_VOL)
                {
                    Speciation(pihm->chemtbl, &pihm->rttbl, 0,
                        &CD->Vcele[i].chms);
                }
            }
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < CD->NumVol; i++)
    {
        int             j;

        for (j = 0; j < pihm->rttbl.NumStc; j++)
        {
            CD->Vcele[i].chms.log10_pconc[j] = log10(CD->Vcele[i].chms.p_conc[j]);
        }
        for (j = 0; j < pihm->rttbl.NumSsc; j++)
        {
            CD->Vcele[i].chms.log10_sconc[j] = log10(CD->Vcele[i].chms.s_conc[j]);
        }
    }

    /* Flux for RIVER flow */
    for (i = 0; i < nriver; i++)
    {
        if (pihm->river[i].down < 0)
        {
            CD->riv = pihm->river[i].wf.rivflow[1];
        }
    }

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        CD->rivd = CD->riv / (DAYINSEC / pihm->ctrl.stepsize);
        //CD->riv = 0;
    }

    for (k = 0; k < pihm->rttbl.NumBTC; k++)
    {
        int             j;

        for (j = 0; j < pihm->rttbl.NumStc; j++)
        {
            if (pihm->rttbl.BTC_loc[k] < 0)
            {
                if (-pihm->rttbl.BTC_loc[k] >= -pihm->rttbl.pumps[0].Pump_Location &&
                    j == pihm->rttbl.pumps[0].Position_Species)
                {
                    CD->Vcele[2 * nelem -pihm->rttbl.BTC_loc[k] - 1].chms.btcv_pconc[j] =
                        log10((CD->Vcele[2 * nelem -pihm->rttbl.BTC_loc[k] - 1].chms.p_conc[j] * CD->riv +
                        pihm->rttbl.pumps[0].Injection_conc * 1.0E-3 * pihm->rttbl.pumps[0].flow_rate) /
                        (CD->riv + pihm->rttbl.pumps[0].flow_rate));
                }
                else
                {
                    CD->Vcele[2 * nelem -pihm->rttbl.BTC_loc[k] - 1].chms.btcv_pconc[j] =
                        CD->Vcele[2 * nelem -pihm->rttbl.BTC_loc[k] - 1].chms.log10_pconc[j];
                }
            }
            else
            {
                CD->Vcele[pihm->rttbl.BTC_loc[k] - 1].chms.btcv_pconc[j] =
                    CD->Vcele[pihm->rttbl.BTC_loc[k] - 1].chms.log10_pconc[j];
            }
        }
    }
}
