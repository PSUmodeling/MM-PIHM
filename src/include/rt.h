typedef struct rtic_struct
{
    double          t_conc[MAXSPS];
    double          p_para[MAXSPS];
} rtic_struct;

typedef struct chmstate_struct
{
    double          t_conc[MAXSPS];         /* concentration (M kg-1 water) */
    double          p_conc[MAXSPS];         /* primary concentration
                                             * (M kg-1 water) */
    double          s_conc[MAXSPS];         /* secondary concentration
                                             * (M kg-1 water) */
    double          p_actv[MAXSPS];         /* activity of primary species */
} chmstate_struct;

typedef struct vol_conc_type
{
    int             index;      /* Volume No. Note this number may be different than the element No. in PIHM */
    int             illness;    /* black list */
    int             type;       /* type of volume: unsaturated, groundwater, river, or river bed */
    double         *t_mole;
    double         *transp_flux;
    double         *react_flux;
    double         *p_para;     /* parameters of primary species
                                 * for surface complexation
                                 * for cation exchage
                                 * for minerals
                                 */
    double          sat;        /* Saturation of Volume, dimensionless, from 0 to 1 */
    double          height_o;   /* height of volume, at previous time step */
    double          height_t;   /* height of volume, at current time step */
    double          height_int; /* temporary variable for intrapolation of gw height */
    double          area;       /* area of the triangular element */
    double          vol_o;      /* volume of water of last time step */
    double          vol;        /* volume of water of current time step */
    double          porosity;   /* porosity of the volume */
    double          rt_step;    /* rt_step of cell (s) */
    double         *log10_pconc;     /* for output only */
    double         *log10_sconc;     /* for output only */
    double         *btcv_pconc; /* for btcv output only */
    rtic_struct     ic;
    chmstate_struct chms;
} vol_conc;

typedef struct face_type
{
    int             nodeup;     /* Volume No. of first node */
    int             nodelo;     /* Volume No. of second node */
    int             nodeuu;     /* Volume No. of node before first node */
    int             nodell;     /* Volume No. of node after second node */
    int             node_trib;  // # node of tributary; > 0 indicates a tributary inflow, 01.14 by Wei Zhi
    int             BC;         /* This face is a boundary face, do not do dispersion and diffusion */
    double          distance;   /* distance from centroid of first node to second node, in meter */
    double          velocity;   /* linear velocity of flux at this face (m s-1) */
    double          flux;       /* flux at surface (m3 s-1) */
    double          flux_trib;  // tributary flux (m3 s-1)
    double          s_area;     /* contact surface area, in square m */
} face;

typedef struct Chem_Data_structure
{
    int             NumVol;     /* Number of total volume in the rt simulator */
    int             NumFac;     /* Number of faces in the rt simulator        */
    double          rivd;       /* Stream discharge in cubic meter per day at the moment */
    double          riv;
    vol_conc       *Vcele;      // An array that stores the volumetric (vol) and chemical (conc) information of grid blocks
    face           *Flux;       // connections between grid blocks
} *Chem_Data;
