// An AdH SuperModel
#ifndef H_SMODEL_SUPER_
#define H_SMODEL_SUPER_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    int myid, npes; // multi-core ID and ncpus
    SFLAGS flags;   // flags associated with superModel
    int id;         // superModel ID
    int itrns;      // the transport constituent ID (CJT -- DO WE NEED THIS?)
    int isSimple;   // Does it only one physics type on whole grid
    char *type;     // surface water, groundwater, transport
    SGRID *grid;    // just a pointer to the grid in the design model
    SDT *ts;        // stores the time information of the model
    int nsubsteps;  // integer divsion of total dt
    int meshcode;   // entire mesh (0), surface(1), or floor(2)
    int forward_step;
    char filebase[MAXLINE]; // base filename for superModel input files

    // superModel solver 
    SLIN_SYS *lin_sys; //pointer to the design model's linear system
    double perturbation;
    double inc_nonlin;
    double tol_nonlin;
    int max_nonlin_linesearch_cuts;
    int max_nonlin_it;
    int it_count_nonlin;
    int force_nonlin_it;
    int nonlinear_it_total;
    int LINEAR_PROBLEM;
    int it_count_nonlin_failed;
    int nalloc_inc;

    // SuperModel Degrees of Freedom
    // Local to process, this should be same as local_range[1]-local_range[0]
    int *my_ndofs;      //pointers to design model, not arrays
    int *my_ndofs_old;
    int *ndofs;         // local number of degrees of freedom
    int *ndofs_old;     //local numer of solution variables the processor is in charge of
    int *macro_ndofs;
    int *macro_ndofs_old;

    // linear systems solution arrays [ndofs]
    double *sol;
    double *sol_old;
    double *sol_older;

    // physics modules
    SSW  *sw;           // surface water storage
    //SCON *con;          // transport constituent storage
    //SSED *sed;          // sediment constituent storage
    //SGW  *gw;           // groundwater storage

    // cjt :: really, this only needs to be done for the body elements
    int *elem1d_physics_mat;    // [nelems1d] - gives the physics material ID
    int *elem2d_physics_mat;    // [nelems2d] - gives the physics material ID
    int *elem3d_physics_mat;    // [nelems3d] - gives the physics material ID
    int nmat_physics; // # of physics materials 
    SMAT_PHYSICS *mat_physics_elem;    // the elemental physics materials
    SMAT_PHYSICS **mat_physics_node;   // the nodal array of pointers to elem material with highest dof connected

    // This is to know all the indepedent variables
    // solved for in this superModel. The actual vars solved 
    // on each node are stored in mat_node
    // If the var position is UNSET_INT, then the variable is not
    // solved for in this superModel
    SIVAR_POSITION ivar_pos;
    
    // This returns position in sol[ndof]
    // If the return is UNSET_INT, then the variable is not solve for
    // on that node
    int **ivars; //[variable_flag (PERTURB_H,...)][nnodes]

    // Series
    int nseries;
    SSERIES *series_head,        *series_curr;
    SSERIES *series_wind_head,   *series_wind_curr;
    SSERIES *series_wave_head,   *series_wave_curr;
    SSERIES *series_gw_psk_head, *series_gw_psk_curr;   

    int *bc_mask;
    double *dirichlet_data;

    // General variables
    double gravity;
    double density;
    int o_flag;
    double initial_grid_mass;
    double grid_mass_error;    

} SMODEL_SUPER;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void smodel_super_alloc_init(SMODEL_SUPER **smod, int nSuperModels);
void smodel_super_free_array(SMODEL_SUPER *sm, int nSuper);
void smodel_super_free(SMODEL_SUPER *smod);
void smodel_super_read(SMODEL_SUPER *smod);
void smodel_super_read_init(SMODEL_SUPER *sm, char *filebase);
void smodel_super_printScreen(SMODEL_SUPER *smod);

//Mark added for testing
void smodel_super_no_read_simple(SMODEL_SUPER *sm, double* dt_in, double* t_init, double* t_prev,
    double* t_final, int nmat_physics, char elemVarCode[4], int isSimple, SGRID *grid, SLIN_SYS *sys);
//from fe_newton.c
//int fe_newton(struct SMODEL_SUPER* sm);

//TODOD, write this up!
int smodel_super_forward_step(SMODEL_SUPER* sm, int (*ts_fnctn)(SMODEL_SUPER*));
int smodel_super_resid(SMODEL_SUPER* sm, double *rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG, int (*fe_resid)(SMODEL_SUPER *, double *, int, double, int, int, int, int));
void smodel_super_update_dirichlet_data(SMODEL_SUPER *sm);
void smodel_super_prep_sol(SMODEL_SUPER *sm);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif

