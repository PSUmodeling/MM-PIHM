#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

/*
 * Function Declarations
 */
void pihm (char *project, int verbose, int debug, char *output_dir, int first_cycle);
void            initialize (char *, Model_Data, Control_Data , N_Vector);
void            initialize_output (char *, Model_Data, Control_Data , char *);
int             f (realtype, N_Vector, N_Vector, void *);
void            read_alloc (char *, Model_Data, Control_Data );
void ReadRiv (char *simulation, Model_Data DS, Control_Data CS);
void ReadMesh (char *simulation, Model_Data DS, Control_Data CS);
void ReadAtt (char *simulation, Model_Data DS, Control_Data CS);
void ReadSoil (char *simulation, Model_Data DS, Control_Data CS);
void ReadGeol (char *simulation, Model_Data DS, Control_Data CS);
void ReadLC (char *simulation, Model_Data DS, Control_Data CS);
void ReadForc (char *simulation, Model_Data DS, Control_Data CS);
void ReadIbc (char *simulation, Model_Data DS, Control_Data CS);
void ReadPara (char *simulation, Model_Data DS, Control_Data CS);
void ReadCalib (char *simulation, Model_Data DS, Control_Data CS);

int Readable (char *token);
int FindLine (FILE *fid, char *token);
void NextLine (FILE *fid, char *cmdstr);
int CountLine (FILE *fid, int num_arg, ...);
int CountOccurance (FILE *fid, char *token);


void            f_update (realtype, realtype *, void *);    /* YS */
void            FreeData (Model_Data, Control_Data );
void            summary (Model_Data, N_Vector, realtype, realtype); /* YS */
realtype        returnVal (realtype, realtype, realtype, realtype);
realtype        CS_AreaOrPerem (int, realtype, realtype, realtype);
void            OverlandFlow (realtype **, int, int, realtype, realtype, realtype, realtype, realtype);
void            OLFeleToriv (realtype, realtype, realtype, realtype, realtype, realtype **, int, int, realtype);
realtype        avgY (realtype, realtype, realtype);
realtype        effKV (realtype, realtype, realtype, realtype, realtype);
realtype        effKV_new (realtype ksatFunc, realtype elemSatn, int status, realtype macKV, realtype KV, realtype areaF);
realtype        effKH (int, realtype, realtype, realtype, realtype, realtype, realtype);
realtype        FieldCapacity (realtype, realtype, realtype, realtype, realtype);
void            is_sm_et (realtype, realtype, void *, N_Vector);
void            PrintInit (Model_Data, char *);

int macpore_status (realtype ksatFunc, realtype elemSatn, realtype gradY, realtype macKV, realtype KV, realtype areaF);

#endif
