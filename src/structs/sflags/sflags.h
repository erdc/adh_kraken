#ifndef H_SFLAGS_
#define H_SFLAGS_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {
    
    /* Model Flags */
    int *model;    // which models are in this superModel
    
    int *params;  // which params are in this superModel

    int *dvars; // which dvars are in this superModel
    /* moving grid */
    bool MG;
    
    /* baroclinic :: flag to signify what variables affect density */
    /* baroclinic code = 0  :: no density effects */
    /* baroclinic code = 1  :: salinity only affects density */
    /* baroclinic code = 10 :: temperature only affects density */
    /* baroclinic code = 11 :: salinity and temperature affect density */
    bool BAROCLINIC;
    
    int EOS;
    /* 0 Linearlized EOS */
    /* 1 Full Equation EOS */
    
    /* mesh adaption */
    bool GRID_ADAPTION;        /* flag that grid *can* be adapted during the simulation */
    bool ADAPTED_THE_GRID;     /* flag that grid was adapted during a given time-step */
    bool GRID_REFINED;         /* flat that grid was refined */
    bool GRID_UNREFINED;       /* flag the grid was unrefined */
    
    /* time adaption */
    bool TIME_ADAPT;           /* was t_adpt_flag */
    bool TIME_ADAPT_FAIL;      /* was t_fail_flag */
    
    /* transport */
    bool TRANSPORT;
    bool NS2_TRANSPORT;
    bool NS3_TRANSPORT;
    bool SW2_TRANSPORT;
    bool SW3_TRANSPORT;
    bool VORTICITY;
    
    /* sediment */
    bool SEDIMENT;
    bool SEDLIB;
    
#ifdef _ADH_GROUNDWATER
    /* groundwater */
    bool GW_SALINITY;
    bool GW_TRANSPORT;
    bool GW_REACTION;
    bool RAY_TRACING;
    bool SOCKETS;
    bool METEOROLOGY;
#endif
    /* steady state flag */
    bool STEADY_STATE;
    
    /* external libraries */
    bool ICM;
    bool NSM;
    
    /* surface stresses */
    bool WAVE_RADS;
    bool WAVE_STRESS;
    bool WIND_STRESS;
    bool WAVE_STATION;
    bool WIND_STATION;
    
    /* coriolis force */
    bool CORIOLIS;
    
    /* manning's unit flag */
    bool MUC;
    
    /* units (MKE or FPS) */
    bool UNITS;
    
    /* conveyance ?? */
    bool CONVEYANCE;
    
    /* flag whether or not to prbooled the adapted mesh */
    bool PRN_ADPT;
    
    /* flag whether to write output or not */
    bool OUTPUT;
    
    /* flag for solver on initial automatic time determination */
    bool SOLVE_ATF;
    
    /* if TRUE, singular matrix occurred, but don't exit */
    bool UMFPACK_FAIL;
    
    /* ice */
    bool ICE;              /* indicates that ice is included (SW) */
    bool INS;              /* indicates that ice is designated by circular string */
    bool nice_coords;      /* number of ice coordinates */
    
    /* tides */
    bool TIDE;
    
    /* unrefine this ts */
    bool UNREFINE;
    
    /* CStorm WSID flag */
    bool CSTORM_WSID;
    
#ifdef _HDF5
    /* Parallel XDMF output */
    bool PC_FILE_XDMF;
#endif
    
    bool FLUX_WEIGHTED_NORMALS; /* Gajanan gkc adding. ON/OFF. Default OFF. Triggered using "NB FLXNML" in bc file. */
    
} SFLAGS;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void sflags_init(SFLAGS *f);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif

