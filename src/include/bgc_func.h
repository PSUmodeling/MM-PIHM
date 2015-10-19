#ifndef BGC_FUNC_H
#define BGC_FUNC_H

void            BgcRead (char *filename, bgc_struct bgc, pihm_struct pihm);
void            ReadAnnFile (ts_struct *ts, char *fn);
void            ReadBinFile (ts_struct *ts, char *fn, int numele);
void            BgcInit (char *simulation, pihm_struct pihm, lsm_struct noah,
    bgc_struct bgc);
void            MapBgcOutput (char *simulation, bgc_struct bgc, int numele,
    char *outputdir);
void            BgcSpinup (char *simulation, bgc_struct bgc, pihm_struct pihm,
    lsm_struct noah);
void            BgcCoupling (int t, int start_time, pihm_struct pihm,
    lsm_struct noah, bgc_struct bgc);
void            Bgc2Noah (int t, pihm_struct pihm, lsm_struct noah, bgc_struct bgc);
void            DailyBgc (bgc_struct bgc, int numele, int numriv, int t, int simstart,
    const double *naddfrac, int first_balance);

void            daymet (const metarr_struct * metarr, metvar_struct * metv,
    int metday);
void            metarr_init (bgc_struct bgc, pihm_struct pihm, int start_time, int end_time);
void            pihm2metarr (bgc_struct bgc, pihm_struct pihm, int t, int start_time, int end_time);

void            presim_state_init (wstate_struct * ws, cstate_struct * cs,
    nstate_struct * ns, cinit_struct * cinit);
void            make_zero_flux_struct (wflux_struct * wf, cflux_struct * cf,
    nflux_struct * nf);
void            restart_input (wstate_struct * ws, cstate_struct * cs,
    nstate_struct * ns, epvar_struct * epv, restart_data_struct * restart);
void            restart_output (wstate_struct * ws, cstate_struct * cs,
    nstate_struct * ns, epvar_struct * epv, restart_data_struct * restart);
void            firstday (const epconst_struct * epc,
    const cinit_struct * cinit, epvar_struct * epv, cstate_struct * cs,
    nstate_struct * ns);
void            zero_srcsnk (cstate_struct * cs, nstate_struct * ns,
    wstate_struct * ws, summary_struct * summary);
void            precision_control (wstate_struct * ws, cstate_struct * cs,
    nstate_struct * ns);

void            radtrans (const cstate_struct * cs,
    const epconst_struct * epc, metvar_struct * metv, epvar_struct * epv,
    double albedo);
void            maint_resp (const cstate_struct * cs,
    const nstate_struct * ns, const epconst_struct * epc,
    const metvar_struct * metv, cflux_struct * cf, epvar_struct * epv);
void            phenology (const epconst_struct * epc,
    const metvar_struct * metv, phenology_struct * phen, epvar_struct * epv,
    cstate_struct * cs, cflux_struct * cf, nstate_struct * ns,
    nflux_struct * nf);
void            leaf_litfall (const epconst_struct * epc, double litfallc,
    cflux_struct * cf, nflux_struct * nf);
void            froot_litfall (const epconst_struct * epc, double litfallc,
    cflux_struct * cf, nflux_struct * nf);
void            soilpsi (const siteconst_struct * sitec, double soilw,
    double *psi);
void            canopy_et (const metvar_struct * metv,
    const epconst_struct * epc, epvar_struct * epv, wflux_struct * wf);
void            total_photosynthesis (const metvar_struct * metv,
    const epconst_struct * epc, epvar_struct * epv, cflux_struct * cf,
    psn_struct * psn_sun, psn_struct * psn_shade);
void            photosynthesis (psn_struct * psn);
void            decomp (double tsoil, const epconst_struct * epc,
    epvar_struct * epv, cstate_struct * cs, cflux_struct * cf,
    nstate_struct * ns, nflux_struct * nf, ntemp_struct * nt);
void            daily_allocation (cflux_struct * cf, cstate_struct * cs,
    nflux_struct * nf, nstate_struct * ns, epconst_struct * epc,
    epvar_struct * epv, ntemp_struct * nt, const double naddfrac,
    const int spinup);
void            annual_rates (const epconst_struct * epc, epvar_struct * epv);
void            growth_resp (epconst_struct * epc, cflux_struct * cf);
void            daily_carbon_state_update (cflux_struct * cf,
    cstate_struct * cs, int alloc, int woody, int evergreen);
void            daily_nitrogen_state_update (nflux_struct * nf,
    nstate_struct * ns, int alloc, int woody, int evergreen);
void            mortality (const epconst_struct * epc, cstate_struct * cs,
    cflux_struct * cf, nstate_struct * ns, nflux_struct * nf);
void            check_carbon_balance (cstate_struct * cs, double *old_balance,
    int first_balance);
void            check_nitrogen_balance (nstate_struct * ns,
    double *old_balance, int first_balance);
void            csummary (cflux_struct * cf, cstate_struct * cs,
    summary_struct * summary);
double          GetCO2 (ts_struct co2_ts, int t);
double          GetNdep (ts_struct ndep_ts, int t);
void nleaching (bgc_grid *grid, int numele, bgc_river *riv, int numriv);

#endif
