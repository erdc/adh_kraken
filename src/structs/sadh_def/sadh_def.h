#ifndef H_SADH_DEF_
#define H_SADH_DEF_

#define MAX_IVAR_DIM 4
#define MAX_TRNS 15 // per dimension


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
    int dim;
    char name[MAXLINE];
    char subname[MAXLINE];
} SVAR;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {
    
    // ++++++++++++++++++++++++++++++++++
    // ADH INDEPENDENT VARIABLES
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
    // ADH INDEPENDENT VARIABLES
    // ++++++++++++++++++++++++++++++++++
    int n_dvars;
    SVAR *dvar;
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
    int _BED_ELE;
    int _GS;
    int _GS_OLD;
    int _DPRS;
    int _DPRS_PLUS;
    int _DPRS_MINUS;
    int _VNODE_FLUX;
    int _WDFLAG;
    int _ERROR;
    
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
    // PHYSICS RESIDUALS
    // ++++++++++++++++++++++++++++++++++
    
    
    
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
