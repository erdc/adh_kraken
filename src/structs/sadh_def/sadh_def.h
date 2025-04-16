#ifndef H_SADH_DEF_
#define H_SADH_DEF_

#define MAX_IVAR_DIM 4
#define MAX_TRNS 15 // per dimension
#define MAX_IVARS 10 //should be updated if ivars are added to sadh_def.c
#define MAX_NDVARS 35 //should be updated if dvars are added to sadh_def.c 
#define MAX_NPARAMS 15
#define ADH_INT 0
#define ADH_DBL 1
#define ADH_SVECT 2
#define ADH_DBL_PTR 3
#define ADH_ELEMENTAL 0
#define ADH_NODAL 1


typedef int (*residual_ptr)(SMODEL_SUPER *, double *, int, double, int, int, int, int);
typedef int (*init_func_ptr)(SMODEL_SUPER *);
residual_ptr select_resid_func(char *model, int *imod, int trns_id);
residual_ptr select_bresid_func(char *model, char *bc_type, char *bc_var, int *imod, int trns_id);
init_func_ptr select_init_func(char *model, int trns_id);


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    int id;
    int nivars;
    char name[MAXLINE];
    char subname[MAXLINE];
    int  var[MAX_IVAR_DIM];
} SADH_MODEL_DEF;

typedef struct {
    int data_type; //Maybe use macros for double, svect2d, int, int * , etc.
    int layout; //For dvar or grid dependent data, elemental vs nodal
    int dim; //spatial dimension of data
    char name[MAXLINE];
    char subname[MAXLINE];
} SVAR;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {
    
    // ++++++++++++++++++++++++++++++++++
    // ADH SOLUTION VARIABLES
    // ++++++++++++++++++++++++++++++++++
    int n_ivars;
    SVAR *ivar;
    int _H;
    int _U;
    int _V;
    int _W;
    int _UDA;
    int _VDA;
    int _DPL;
    int _PRS;
    int _HEAT;
    int _SAL;
    int _CON[MAX_TRNS];
    
    // ++++++++++++++++++++++++++++++++++
    // ADH DEPENDENT VARIABLES
    // grid varying quantities
    // can be elemental or nodal based
    // ++++++++++++++++++++++++++++++++++
    int n_dvars;
    SVAR *dvar;
    //Maybe wind is just 1 and we can 
    //end up allocating this as nodal/elemental pointer?
    int _WIND_SX;
    int _WIND_SY;
    int _WAVE_SX;
    int _WAVE_SY;
    int _WAVE_XX;
    int _WAVE_XY;
    int _WAVE_YX;
    int _WAVE_YY;
    int _DENSITY;
    int _BED_DPL;
    int _BED_DPL_OLD;
    int _BED_ELE;
    int _GS;
    int _GS_OLD;
    int _DPRS;
    int _DPRS_PLUS;
    int _DPRS_MINUS;
    int _DPL_PERTURBATION;
    int _CONT_ERROR;
    int _HYD_VISCOSITY;
    int _TRN;
    int _GRAD_BED_X;
    int _GRAD_BED_Y; 
    //tangent vector
    int _TAN_VEC_X;
    int _TAN_VEC_Y;
    int _SURFACE_VEL_X; // SVECT2D surface_vel;
    int _SURFACE_VEL_Y; // 
    int _BOTTOM_VEL_X; //
    int _BOTTOM_VEL_Y; //    //SVECT2D bottom_vel;

    int _VNODE_FLUX;
    int _WDFLAG;
    int _ERROR;
    int _ELEM_RHS_SUPG_DACONT; //for SW3 only
    int _ELEM_RHS_SUPG_CONT; // for SW3 only
    int _ELEM_RHS_DACONT_EXTRA_TERMS; // for SW2 only, FLIPPING ORDER FOR EFFICIENCY [nsw_elems][nnode_on_elem]

    
    // ++++++++++++++++++++++++++++++++++
    // ADH MODELS
    // ++++++++++++++++++++++++++++++++++
    int nmodels;
    SADH_MODEL_DEF *model;
    int _SW1, _SW2, _SW3;
    int _NS3;
    int _NS3SPLIT;
    int _DW2D;
    int _WVEL;
    int _PRESSURE;
    int _POISSON;
    int _HEAT2D;
    int _GW3D;
    int _TRNS1D[MAX_TRNS], _TRNS2D[MAX_TRNS], _TRNS3D[MAX_TRNS];
    
    
    // ++++++++++++++++++++++++++++++++++
    // ADH PARAMETERS
    // Any possible constants on grid
    // that could be declared
    // and set in a bc file
    // Parmeters are doubles or ints (use type)
    // ++++++++++++++++++++++++++++++++++
    int n_params;
    SVAR *param;
    //Can have multiple instances such as SSW_params etc
    //but for now we will try to keep it all in one
    int _WD_L_LIMIT;                                                 
    int _WD_U_LIMIT;
    int _WD_L_TOL;
    int _WD_U_TOL;
    int _WD_RATE_L_TOL;
    int _WD_RATE_U_TOL;
    int _VISCOSITY;
    int _MANNINGS_UNITS_CONSTANT;
    int _TAU_PG;
    int _GRAVITY;
    int _DENSITY_PARAM;
    int _O_FLAG;
    int _INIT_GRID_MASS;
    int _GRID_MASS_ERROR;   
    int _ELEM_RHS_REALLOC; //underlying data is integer

    
    
} SADH_DEF;

extern SADH_DEF adh_def;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// Methods
void sadh_def_init(void);
void sadh_def_sivar_printScreen(SVAR v);
void sadh_def_model_printScreen(SADH_MODEL_DEF *m);

//int (*select_resid_func(char *model, int *imod))(struct SMODEL_SUPER *, double *, int, double, int, int, int, int);
//int (*select_bresid_func(char *model, char *bc_type, char *bc_var, int *imod))(struct SMODEL_SUPER *, double *, int, double, int, int, int, int);
//int (*select_init_func(char *model))(struct SMODEL_SUPER *);
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



#endif
