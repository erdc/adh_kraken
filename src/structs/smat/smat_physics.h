#ifndef H_SMAT_PHYSICS_
#define H_SMAT_PHYSICS_

typedef struct SMODEL_SUPER SMODEL_SUPER;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    int imod;    // the SADH_DEF integer ID for model data associated with pde
    int iseries; // hold time-series ID when needed
    int (*resid)(SMODEL_SUPER *, double *, int, double, int, int, int, int);
    int (*init)(SMODEL_SUPER *); // need to point this
} SPDE;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    int id;         // Model ID
    int n;          // number of independent variables on this material
    int ntrns;      // # of transport constituents on this material
    SIVAR_POSITION ivar_pos; // a list of independent variables on the element
    int *ivar_loc;  // array of n who's entries are location in **ivar array
    
    int npdes;      // # of partial differential equations/residuals on this mat
    SPDE *pde;      // stores info on the PDEs/residuals
} SMAT_PHYSICS;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// Methods
void smat_physics_alloc_init(SMAT_PHYSICS *mat);
void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat);
void smat_physics_free_array(SMAT_PHYSICS *mat, int nmat);
void smat_physics_position_flag(SMAT_PHYSICS **mat_node, int nnodes, int *FLAG);
void smat_physics_printScreen(SMAT_PHYSICS *m);
void smat_physics_update_array(SMAT_PHYSICS *m, int nmat_physics, SIVAR_POSITION* ivar_pos);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
