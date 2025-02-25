// An AdH DESIGN MODEL
#ifndef H_SMODEL_DESIGN_
#define H_SMODEL_DESIGN_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    int model1;
    int model2;
    int type;
} SMODEL_COUPLE;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    char filebase[MAXLINE];        // the base project filename
    char filename_design[MAXLINE]; // file - design data
    char filename_grid[MAXLINE];   // file - design grid
    char filename_params[MAXLINE]; // file - supermodel parameters
    FILE *xmf;

    int nSuperModels;           // total # of superModels
    SMODEL_SUPER *superModel;   // superModel array
    SMODEL_COUPLE *superCouple; // superModel coupling type
    int nMono;   // number of monolothic-coupled supermodels
    int nSimple; // number of single physics supermodels (nMono+nSimple = nSuperModels)

    SGRID *grid;  // the designer grid
    SDT ts;       // top-level time-stepping info

    SSERIES *series_dt;  // time-step series
    SSERIES *series_out; // output series

    SCOVERAGE *params;   // superModel shared parameter coverage


    // only allocate memory for linear systems of each sparsity pattern
    int nlin_sys; // number of unique supermodels, leading to unique sparsity structure
                    // a Unique supermodel is one that is either monolithically coupled or
                    // unique sparsity
    int *lin_sys_id; // array of [nSuperModels] that gives the index of the linear system it belongs to
                     //likely unnecessary   
    int *unique_id; // integer array of size nlin_sys, saves index of first Unique super model for each sparsity pattern

    SLIN_SYS *lin_sys; // array of [nlin_sys] systems

    // int *nFluxInterfaces; // total # of flux interfaces between supermodesl
    // SINTERFACE_FLUX  will hold all data for interfacing between supermodels
    
} SMODEL_DESIGN;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void smodel_design_defaults(SMODEL_DESIGN *dm);
void smodel_design_alloc_init(SMODEL_DESIGN **dmod, int nSuperModels);
void smodel_design_free(SMODEL_DESIGN *dmod);
void smodel_design_printScreen(SMODEL_DESIGN *dmod);
#ifdef _MESSG
int smodel_design_init(SMODEL_DESIGN *dmod, char *filename, bool input_check, MPI_Comm *comm);
#else
int smodel_design_init(SMODEL_DESIGN *dmod, char *filename, bool input_check);
#endif
void smodel_design_read(SMODEL_DESIGN *dmod, char *filename);
void smodel_design_init_no_read(SMODEL_DESIGN *dmod, double dt_in, double t_init, double t_final,
    int nSuperModels, int nphysics_mats[], char ***elemVarCode, int **coverage_arrays);
void smodel_design_xmf_init(SMODEL_DESIGN *dm, char *filename, char *domain_name);
void smodel_design_xmf_write(SMODEL_DESIGN *dm, int mesh_no);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
