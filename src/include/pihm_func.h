#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

/*
 * Function Declarations
 */
void PIHMRun (char *simulation, char *outputdir, int first_cycle);
void CreateOutputDir (char *project, char *outputdir, int overwrite_mode);
//void            initialize (char *, Model_Data, Control_Data , N_Vector);
//void            initialize_output (char *, Model_Data, Control_Data , char *);
//int             f (realtype, N_Vector, N_Vector, void *);
void ReadAlloc (char *simulation, pihm_struct pihm);
void ReadRiv (char *project, riv_att_tbl_struct *riv_att_tbl, riv_shp_tbl_struct *riv_shp_tbl, riv_matl_tbl_struct *riv_matl_tbl, riv_ic_tbl_struct *riv_ic_tbl, ic_struct *ic, forcing_ts_struct *forcing);
void ReadMesh (char *project, mesh_tbl_struct *mesh_tbl);
void ReadAtt (char *project, attrib_tbl_struct *attrib_tbl, ic_struct *ic, int numele);
void ReadSoil (char *project, soil_tbl_struct *soil_tbl);
void ReadGeol (char *project, geol_tbl_struct *geol_tbl);
void ReadLC (lc_tbl_struct *lc_tbl);
void ReadForc (char *project, forcing_ts_struct *forcing);
void ReadLAI (char *project, forcing_ts_struct *forcing, int numele, const attrib_tbl_struct *attrib_tbl);
void ReadIbc (char *project, forcing_ts_struct *forcing);
void ReadPara (char *project, ctrl_struct *ctrl);
void ReadCalib (char *project, char *simulation, calib_struct *cal);
void ReadInit (char *project, char *simulation, ic_struct *ic, int numele, int numriv);
void            FreeData (pihm_struct pihm);

//void ReadCalib (char *simulation, Model_Data DS, Control_Data CS);
//
int Readable (char *token);
int FindLine (FILE *fid, char *token);
void NextLine (FILE *fid, char *cmdstr);
int CountLine (FILE *fid, int num_arg, ...);
void ReadTS (char *cmdstr, int *ftime, double *data, int nvrbl);
void CheckFile (FILE *fid, char *fn);
void ReadKeywordDouble (char *buffer, char *keyword, double *value);
void ReadKeywordInt (char *buffer, char *keyword, int *value);
void ReadKeywordTime (char *buffer, char *keyword, int *value);
//int CountOccurance (FILE *fid, char *token);
//
//
void Initialize (pihm_struct pihm, N_Vector CV_Y);
void InitMeshStruct (elem_struct *elem, int numele, mesh_tbl_struct mesh_tbl);
void InitTopo (elem_struct *elem, int numele, mesh_tbl_struct mesh_tbl);
void InitSoil (elem_struct *elem, int numele, attrib_tbl_struct attrib_tbl, soil_tbl_struct soil_tbl, geol_tbl_struct geol_tbl, calib_struct cal);
double FieldCapacity (double alpha, double beta, double kv, double thetas, double thetar);
double WiltingPoint (double thetas, double thetar, double alpha, double beta);
void InitLC (elem_struct *elem, int numele, attrib_tbl_struct attrib_tbl, lc_tbl_struct lc_tbl, calib_struct calibration);
void InitRiver (river_struct *riv, int numriv, elem_struct *elem, riv_att_tbl_struct riv_att_tbl, riv_shp_tbl_struct riv_shp_tbl, riv_matl_tbl_struct riv_matl_tbl, mesh_tbl_struct mesh_tbl, calib_struct cal);
void InitForcing (elem_struct *ele, int numele, attrib_tbl_struct attrib_tbl, forcing_ts_struct *forcing, calib_struct cal);
void CorrectElevation (elem_struct *elem, int numele, river_struct *riv, int numriv);
void InitSurfL (elem_struct *ele, int numele, river_struct *riv, mesh_tbl_struct mesh_tbl);
void SaturationIC (const elem_struct *ele, int numele, const river_struct *riv, int numriv, ic_struct *ic);
void InitStateVrbl (elem_struct *elem, int numele, river_struct *riv, int numriv, N_Vector CV_Y, ic_struct ic);
void CalcModelStep (ctrl_struct *ctrl);

void MapOutput (char *simulation, pihm_struct pihm, char *outputdir);
void InitOutputFile (prtctrl_struct *prtctrl, int nprint, int ascii);

//void            f_update (realtype, realtype *, void *);    /* YS */
//void            summary (Model_Data, N_Vector, realtype, realtype); /* YS */
//realtype        returnVal (realtype, realtype, realtype, realtype);
//realtype        CS_AreaOrPerem (int, realtype, realtype, realtype);
//void            OverlandFlow (realtype **, int, int, realtype, realtype, realtype, realtype, realtype);
//void            OLFeleToriv (realtype, realtype, realtype, realtype, realtype, realtype **, int, int, realtype);
//realtype        avgY (realtype, realtype, realtype);
//realtype        effKV (realtype, realtype, realtype, realtype, realtype);
//realtype        effKV_new (realtype ksatFunc, realtype elemSatn, int status, realtype macKV, realtype KV, realtype areaF);
//realtype        effKH (int, realtype, realtype, realtype, realtype, realtype, realtype);
//realtype        FieldCapacity (realtype, realtype, realtype, realtype, realtype);
//void            is_sm_et (realtype, realtype, void *, N_Vector);
//void            PrintInit (Model_Data, char *);
//
//int macpore_status (realtype ksatFunc, realtype elemSatn, realtype gradY, realtype macKV, realtype KV, realtype areaF);
double          KrFunc (double alpha, double beta, double satn);

#endif
