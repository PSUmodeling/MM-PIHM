/**********************************************************************************
 * File        : pihm.h                                                           *
 * Function    : Declaration and Definition of global variables and data structure*
 * Developer of PIHM 2.0: Mukesh Kumar (muk139@psu.edu)                           *
 * Developer of PIHM 1.0: Yizhong Qu   (quyizhong@gmail.com)                      *
 * Version     : Nov, 2007 (2.0)                                                  *
 *--------------------------------------------------------------------------------*
 *                                                                                *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...............................* 
 * a) Definition of new variables for ELEMENT data model related to 		  *
 *	i. Subsurface: KsatH, KsatV, infKsatV,infD, RzD, macD, macKsatH, macKsatV *
 *	vAreaF, hAreaF								  * 
 *	ii. ET: LAImax, vegFrac, Albedo, Rs_ref, Rmin				  *
 *	iii. Snow: meltF							  *
 *	iv. Surface: dhBydx, dhBYdy, surfH					  *
 * b) Definition of new variables for RIVER data model				  *
 *	i. KsatH, KsatV, bedThick 						  *
 * c) Definition of new structures:						  *
 *	i. geology								  *
 *	ii. Land Cover								  *
 *	iii. Calibration							  *
 * d) Definition of New Control Parameters for each file			  *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                	  *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 **********************************************************************************/

//#include "sundialstypes.h"
#include "sundials_types.h"   // 09.16
#include "nvector_serial.h"

/* Definition of Global Variable Types */

FILE *riv_state_file;
typedef struct element_type	/* Data model for a triangular element */ 
	{
  	int index; 		/* Element No. */ 
  	int node[3];        	/* anti-clock-wise */
  	int nabr[3];        	/* neighbor i shares edge i (0: on boundary) */
    
  	realtype edge[3];   	/* edge i is from node i to node i+1 */
  	realtype area;      	/* area of element */
  
  	realtype x;         	/* x of centroid */
  	realtype y;         	/* y of centroid */
  	realtype zmin;      	/* z_min of centroid */
  	realtype zmax;      	/* z_max of centroid */
  
  	realtype KsatH;		/* horizontal geologic saturated hydraulic conductivity */
  	realtype KsatV;		/* vertical geologic saturated hydraulic conductivity */
	realtype infKsatV;	/* vertical surface saturated hydraulic conductivity */
	realtype effKV;		/* Shi */
  	realtype Porosity;
	realtype ThetaS;
	realtype ThetaR;
	realtype ThetaW;
	realtype ThetaRef;		
	realtype infD;		/* depth from ground surface accross which head is calculated during infiltration */
  	realtype Alpha;		/* Alpha from van-genuchten eqn which is given by satn = 1/pow(1+pow(abs(Alpha*psi),Beta),1-1/Beta) */
  	realtype Beta;
  	realtype RzD;		/* Root zone depth */
  	realtype macD;		/* macropore Depth */
  	realtype macKsatH;	/* macropore horizontal saturated hydraulic conductivity */
  	realtype macKsatV;	/* macropore vertical saturated hydraulic conductivity */
  	realtype vAreaF;	/* macropore area fraction on a vertical cross-section */
  	realtype hAreaF;	/* macropore area fraction on a horizontal cross-section */
  	int Macropore;       	/* 1: macropore; 0: regular soil */

	int MPL;		/* The deepest layer that macropore reaches */
	int rootL;		/* The layer that root reaches */

  	realtype LAImax;    	/* maxm. LAI accross all seasons for a vegetation type */
  	realtype VegFrac;   	/* areal vegetation fraction in a triangular element */ 
  	realtype Albedo;    	/* albedo of a triangular element */
	realtype Albedo_min;
	realtype Albedo_max;
	realtype Emiss_min;
	realtype Emiss_max;
	realtype z0_min;
	realtype z0_max;
	realtype h_s;		/* parameter used in canopy resistance calculation */
  	realtype Rs_ref;     	/* reference incoming solar flux for photosynthetically active canopy */
  	realtype Rmin;       	/* minimum canopy resistance */
  	realtype Rough;      	/* surface roughness of an element */
  
  	realtype windH;		/* wind measurement height */
	  realtype tempK; 
  
    
  	int soil;           	/* soil type */
	int geol;		/* geology type */
  	int LC;             	/* Land Cover type  */
  	int IC;             	/* initial condition type */
  	int BC[3];             	/* boundary type. 0:natural bc (no flow); 1:Dirichlet BC; 2:Neumann BC */
  	int prep;           	/* precipitation (forcing) type */
  	int temp;           	/* temperature (forcing) type   */
  	int humidity;       	/* humidity type */
  	int WindVel;        	/* wind velocity type  */
//  	int Rn;             	/* net radiation input */
//  	int G;              	/* radiation into ground */
  	int Sdown;             	/* downward solar radiation, modified by Y. Shi */
  	int Ldown;              /* downward long wave radiation, modified by Y. Shi */
  	int pressure;       	/* pressure type */
  	int source;         	/* source (well) type */
 	int meltF;		/* meltFactor */ 
	/* for calculation of dh/ds */
  	realtype surfH[3];	/* Total head in neighboring cells */	
  	realtype surfX[3];	/* Center X location of neighboring cells */
  	realtype surfY[3];	/* Center Y location of neighboring cells */
  	realtype dhBYdx;	/* Head gradient in x dirn. */
  	realtype dhBYdy;	/* Head gradient in y dirn. */
	} element;

typedef struct nodes_type	/* Data model for a node */
	{
  	int index; 		/* Node no. */
     
  	realtype x;          	/* x coordinate */
  	realtype y;          	/* y coordinate */
  	realtype zmin;       	/* z bed rock elevation */
  	realtype zmax;       	/* z surface elevation  */
  
	} nodes;

typedef struct element_IC_type	/* Initial state variable conditions on each element */
	{
  	int index;
  
  	realtype interception;	/* Interception storage (Note all these variables have dimension of L */
  	realtype snow;		/* Snow depth */
  	realtype surf;		/* Overland flow depth */
  	realtype unsat;		/* unsaturated zone depth */
  	realtype sat;		/* saturated zone depth */
  
	} element_IC;

typedef struct soils_type
	{
  	int index;           	/* index */
  
  	realtype KsatV;       	/* vertical saturated soil conductivity */
  	realtype ThetaS;      	/* soil porosity */
  	realtype ThetaR;      	/* soil moisture residual */
	realtype ThetaW;	/* wilting point, expanded by Y. Shi */
	realtype ThetaRef;	/* field capacity, expanded by Y. Shi */
  	realtype Alpha;      	/* soil curve parameter 1 */
  	realtype Beta;       	/* soil curve parameter 2 */
  
  	realtype hAreaF;       	/* macroporous area fraction on horizontal section */
  	realtype macKsatV;      /* macroporous saturated vertical conductivity */
 
  	realtype infD;		/* depth from ground surface accross which head is calculated during infiltration */
	} soils;

typedef struct geol_type
        {
        int index;              /* index */

        realtype KsatH;         /* horizontal saturated geology conductivity */
        realtype KsatV;         /* vertical saturated geology conductivity */
        realtype ThetaS;        /* geology porosity */
        realtype ThetaR;        /* residual porosity */
        realtype Alpha;         /* van genuchten parameter */
        realtype Beta;          /* van genuchten parameter */

        realtype vAreaF;        /* macroporous area fraction on vertical section */
        realtype macKsatH;      /* macroporous saturated horizontal conductivity */
        realtype macD;

        } geol;

typedef struct lc_type
	{
  	int index;           	/* index */

  	realtype LAImax;     	/* max LAI */
  	realtype VegFrac;    	/* Canopy Fracn */
  	realtype Albedo_min;     	/* Minimum albedo, added by Y. Shi */
  	realtype Albedo_max;     	/* Maximum albedo, added by Y. Shi */
  	realtype Emiss_min;     	/* Minimum emissivity, added by Y. Shi */
  	realtype Emiss_max;     	/* Maximum emissivity, added by Y. Shi */
  	realtype z0_min;     	/* Minimum roughness length, added by Y. Shi */
  	realtype z0_max;     	/* Maximum roughness length, added by Y. Shi */
	realtype h_s;		/* Vapor pressure deficit stress parameter, added by Y. Shi */
  	realtype Rs_ref;     	/* Visible solar flux used in radiation stress */
  	realtype Rmin;       	/* Minimum stomatal resistance */
  	realtype Rough;      	/* Surface roughness factor  */
  	realtype RzD;	       	/* RootZone Depth*/
	} LC;

typedef struct river_segment_type
	{
  	int index;
  
  	realtype x;          	/* x of river segment */
  	realtype y;  
  	realtype zmin;       	/* bed elevation  */
  	realtype zmax;       	/* bank elevation */
  	realtype depth;      	/* max depth */
  	realtype Length;	/* Riv segment Length */
  	realtype Rough;		/* Manning's roughness coeff */
 	realtype KsatH;		/* Side conductivity */
	realtype KsatV;		/* Bed conductivity */
	realtype bedThick;
	realtype coeff;		/* Coefficient c in D = c*pow(B/2,interpOrd) where D is depth*/ 
  	int FromNode;		/* Upstream Node no. */
  	int ToNode;		/* Dnstream Node no. */
  	int down;            	/* down stream segment */
  	int LeftEle;		/* Left neighboring element */
  	int RightEle;		/* Right neighboring element */
  	int shape;           	/* shape type    */
  	int material;        	/* material type */
  	int IC;              	/* IC type */
  	int BC;              	/* BC type */
  	int reservoir;
  
	} river_segment;

typedef struct river_shape_type
	{
  	int index;           
  	realtype depth;      	/* depth */
  	int interpOrd;       	/* Interpolation order for river shape: Order =1 (rectangle), 2(triangle), 3(quadratic) and 4(cubic)*/
  	realtype coeff;	       	/* Coefficient c in D = c*pow(B/2,interpOrd) */
  
	} river_shape;

typedef struct river_material_type
	{
  	int index;
  	realtype Rough;
  	realtype Cwr;		/* Weir Discharge Coefficient */
 	realtype KsatH;		/* Conductivity of river banks */
	realtype KsatV;		/* Conductivity of river bed */
	realtype bedThick; 	/* thickeness of conductive river bed */
	} river_material;

typedef struct river_IC_type
	{
  	int index;
  	realtype value;		/* initial flow depth */
  
	} river_IC;

typedef struct TSD_type
	{
  	char name[5];
  	int index;
  	int length;	/* length of time series */
  	int iCounter;	/* interpolation counter */
  	realtype **TS;	/* 2D time series data */
  
	} TSD;

typedef struct global_calib
	{
        realtype KsatH;		/* For explanation of each calibration variable, look for corresponding variables above */
        realtype KsatV;
        realtype infKsatV;
        realtype macKsatH;
        realtype macKsatV;
        realtype infD;
        realtype RzD;
        realtype macD;
        realtype Porosity;
        realtype Alpha;
        realtype Beta;
        realtype vAreaF;
        realtype hAreaF;
	realtype Temp;
	realtype Prep;
        realtype VegFrac;
        realtype Albedo;
        realtype Rough;

        realtype rivRough;
        realtype rivKsatH;
        realtype rivKsatV;
        realtype rivbedThick;
	realtype rivDepth;
	realtype rivShapeCoeff;

	realtype TF;
	realtype IS;
	realtype Rmin;
	realtype Czil;
	realtype fx_soil;
	realtype fx_canopy;
	realtype Rs_ref;
	realtype h_s;
	realtype Tref;
	realtype ThetaRef;
	realtype ThetaW;

	  realtype PCO2;
	  realtype Keq;
	  realtype SSA;
	  realtype Site_den;
	  realtype Prep_conc;
	  realtype GW_tot;
	  realtype GW_conc;
	} globalCal; 

typedef struct process_control
	{
	realtype Et0;
	realtype Et1;
	realtype Et2;
	} processCal;

typedef struct model_data_structure 	/* Model_data definition */
	{
  	int UnsatMode;               	/* Unsat Mode */
  	int SurfMode;                	/* Surface Overland Flow Mode */
  	int RivMode;                 	/* River Routing Mode */
    
  	int NumEle;                  	/* Number of Elements */
  	int NumNode;                 	/* Number of Nodes    */
  	int NumRiv;                  	/* Number of Rivere Segments */

  	int NumPrep;                 	/* Number of Precipatation time series types  */
  	int NumTemp;                 	/* Number of Temperature time series types      */
  	int NumHumidity;             	/* Number of Humidity time series types         */
  	int NumWindVel;              	/* Number of Wind Velocity time series types    */
//  	int NumRn;                   	/* Number of Net Radiation time series types    */
//  	int NumG;                    	/* Number of Ground Heat time series types      */
  	int NumSdown;                  	/* Number of downward solar radiation time series types, modified by Y. Shi    */
  	int NumLdown;                  	/* Number of downward longwave radiation time series types, modified by Y. Shi      */
  	int NumP;                    	/* Number of Pressure time series types         */
  	int NumSource;               	/* Number of Source time series types           */
  
  	int NumSoil;                 	/* Number of Soils           */
  	int NumGeol;                 	/* Number of Geologies           */
  	int NumRes;                  	/* Number of Reservoir       */
  	int NumLC;                   	/* Number of Land Cover Index Data */

	realtype **dsoil;		/* Soil layer depths for soil temperature calculation, modified by Y. Shi */

  	int NumMeltF;		       	/* Number of Melt Factor Time series */
  
  	int Num1BC;                  	/* Number of Dirichlet BC    */
  	int Num2BC;                  	/* Number of Numann BC       */
  	int NumEleIC;                	/* Number of Element Initial Condtion */
  
  	int NumRivShape;             	/* Number of River Shape     */
  	int NumRivMaterial;          	/* Number of River Bank/Bed Material */
  	int NumRivIC;                	/* Number of River Initial Condition  */
  	int NumRivBC;                	/* Number of River Boundary Condition */
 
	realtype TF;
	realtype IS;
	realtype Czil;
	realtype fx_soil;
	realtype fx_canopy;
	realtype Tref;

  	element *Ele;                	/* Store Element Information  */
  	nodes *Node;                 	/* Store Node Information     */
  	element_IC *Ele_IC;          	/* Store Element Initial Condtion */
  	soils *Soil;                 	/* Store Soil Information     */
  	geol *Geol;                 	/* Store Geology Information     */
  	LC *LandC;		       	/* Store Land Cover Information */
  
  	river_segment *Riv;          	/* Store River Segment Information */
  	river_shape *Riv_Shape;      	/* Store River Shape Information   */
  	river_material *Riv_Mat;     	/* Store River Bank Material Information */
  	river_IC *Riv_IC;            	/* Store River Initial Condition   */

  	TSD *TSD_Inc;                	/* Infiltration Capacity Time Series Data */
  	TSD *TSD_LAI;                	/* Leaves Area Index Time Series Data     */
//  	TSD *TSD_DH;		     	/* Zero plane Displacement Height */  
  	TSD *TSD_RL;		     	/* Roughness Length */  
  	realtype *ISFactor;          	/* ISFactor is used to calculate ISMax from LAI */
	realtype *windH;		/* Height at which wind velocity is measured */
  	TSD  *TSD_MeltF;	       	/* Monthly Varying Melt Factor for Temperature Index model */  

  	TSD *TSD_EleBC;              	/* Element Boundary Condition Time Series Data  */
  	TSD *TSD_Riv;                	/* River Related Time Series Data  */
  	TSD *TSD_Prep;               	/* RainFall Time Series Data       */
  	TSD *TSD_Temp;               	/* Temperature Time Series Data    */
  	TSD *TSD_Humidity;           	/* Humidity Time Series Data       */
  	TSD *TSD_WindVel;            	/* Wind Velocity Time Series Data  */
//  	TSD *TSD_Rn;                 	/* Net Radiation Time Series Data  */
//  	TSD *TSD_G;                  	/* Radiation into Ground Time Series Data */
  	TSD *TSD_Sdown;                	/* Net Radiation Time Series Data, modified by Y. Shi  */
  	TSD *TSD_Ldown;                	/* Radiation into Ground Time Series Data, modified by Y. Shi */
  	TSD *TSD_Pressure;           	/* Vapor Pressure Time Series data       */
  	TSD *TSD_Source;             	/* Source (well) Time Series data  */

  	realtype **FluxSurf;     	/* Overland Flux   */
  	realtype **FluxSub;      	/* Subsurface Flux */
  	realtype **FluxRiv;      	/* River Segement Flux */
	  realtype **VeloSub;             /* velocity of subsurface fluxes modified by C. Bao*/
	  realtype **DistSub;             /* Distance of controids, modified by C. Bao*/
	  realtype **AreaSub;             /* area of the contact surface, modified by C. Bao */

  	realtype *ElePrep;		/* Precep. on each element */
  	realtype *EleETloss;		
  	realtype *EleNetPrep;		/* Net precep. on each elment */
  	realtype *EleViR;		/* Variable infiltration rate */
  	realtype *Recharge;		/* Recharge rate to GW */
  	realtype *EleSnow;		/* Snow depth on each element */
  	realtype *EleSnowGrnd;		/* Snow depth on ground element */
  	realtype *EleSnowCanopy;	/* Snow depth on canopy element */
  	realtype *EleIS;		/* Interception storage */
  	realtype *EleISmax;		/* Maximum interception storage (liquid precep) */
	realtype *EleISsnowmax;		/* Maximum interception storage (snow) */
  	realtype *EleTF;		/* Through Fall */
  	realtype **EleET;		/* Evapo-transpiration (from canopy, ground, subsurface, transpiration) */
	realtype *EleDew;		/* Dew, expanded by Y. Shi */
	realtype *EleEp;		/* Potential ET, expanded by Y. Shi */
	realtype *EleH;			/* Sensible heat flux, expanded by Y. Shi */
	realtype *EleLE;		/* Latent heat flux, expanded by Y. Shi */
	realtype *EleG;			/* Ground heat flux, expanded by Y. Shi */
	realtype *EleOVLbuffer;		/* Stores surface water level of last time step, expanded by Y. Shi */
	realtype *EleGWbuffer; 		/* Stores groundwater level of last time step, expanded by Y. Shi */
	realtype *Eleunsatbuffer; 	/* Stores unsaturated storage of last time step, expanded by Y. Shi */
	realtype *C_h;			/* exchange coefficient of heat, expanded by Y. Shi */
	realtype *C_m;			/* exchange coefficient of mass, expanded by Y. Shi */
	realtype **EleSM;		/* Soil volumetric water content, expanded by Y. Shi */
	realtype *Tsfc;			/* Surface temperature, expanded by Y. Shi */
	realtype *Tbot;			/* Soil temperatrure of the bottom layer, expanded by Y. Shi */
	int *EleBot;		/* Bottom layer of soil moisture model, expanded by Y. Shi */
	realtype **Tsoil;		/* Soil temperatrures, expanded by Y. Shi */
  	realtype Q; 
 	realtype *DummyY; 
//	realtype *PrintVar[22];
	realtype *PrintVar[35];		/* Modified by Y. Shi */
	processCal pcCal;
	} *Model_Data;

typedef struct control_data_structure
	{
  	int Verbose;
  	int Debug;
  
  	int Solver;	/* Solver type */
  	int NumSteps;	/* Number of external time steps (when results can be printed) for the whole simulation */	
  
  	int gwD;	/* File boolean, Choose 1 if your want to print ground water */
	int surfD;	/* File boolean, Choose 1 if your want to print overland flow */
	int snowD;      /* File boolean, Choose 1 if your want to print snow Depth */
	int rivStg;     /* File boolean, Choose 1 if your want to print river Stage */
	int Rech;       /* File boolean, Choose 1 if your want to print recharge to ground water */
	int IsD;        /* File boolean, Choose 1 if your want to print interception depth */
	int usD;        /* File boolean, Choose 1 if your want to print unsaturated depth */
	int lsv;	/* File boolean, Choose 1 if you want to print land surface variables */
	int et[3];      /* File boolean, Choose 1 if your want to print individual evapo-transpiration components */
	int rivFlx[10];	/* File boolean, Choose 1 if your want to print river/river bed fluxes */
  
  	int gwDInt;	/* Time interval to output average val of variables */
	int surfDInt;
	int snowDInt;
	int rivStgInt;
	int RechInt;
	int IsDInt;
	int usDInt;
	int etInt;
	int rivFlxInt;
  
  	int init_type;		/* initialization mode */
  
  	realtype abstol;	/* absolute tolerance */
  	realtype reltol;	/* relative tolerance */
  	realtype InitStep;	/* initial step size */
  	realtype MaxStep;	/* Maximum step size */
  	realtype ETStep;	/* Step for et from interception */
  
  	int GSType, MaxK;	/* Maximum Krylov order */
  	realtype delt;
  
  	realtype StartTime;	/* Start time of simulation */
  	realtype EndTime;	/* End time of simulation */
  
  
  	int outtype;		
  	realtype a;		/* External time stepping controls */
  	realtype b;
  
  	realtype *Tout;  
 	
	globalCal Cal;		/* Convert this to pointer for localized calibration */ 
	} Control_Data;

// 09.23
//void FPrintFinalStats(FILE *, long int iopt[], realtype ropt[]);
//void PrintData(FILE **,Control_Data *, Model_Data, N_Vector, realtype);

