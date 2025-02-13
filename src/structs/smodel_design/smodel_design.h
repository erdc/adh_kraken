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

    char filename_design[MAXLINE]; // file - design data
    char filename_grid[MAXLINE];   // file - design grid
    char filename_params[MAXLINE]; // file - supermodel parameters

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

    // the number of unique sparsity patterns and which superModel has them
    int nUnique;    // number of unique supermodels, leading to unique sparsity structure
                    // a Unique supermodel is one that is either monolithically coupled or
                    // unique sparsity
    int *unique_id; // integer array of size nUnique, saves index of first Unique super model for each sparsity pattern

    // only allocate memory for linear systems of each sparsity pattern
    int nlin_sys;
    int *lin_sys_id; // array of [nSuperModels] that gives the index of the linear system it belongs to
    SLIN_SYS *lin_sys; // array of [nUnique] systems

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
void smodel_design_no_read_simple(SMODEL_DESIGN *dm, double dt_in, double t_init, double t_final, int nphysics_mat, char elemVarCode[4] ,SGRID *grid);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
