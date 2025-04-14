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
    //int n;          // number of independent variables on this material, in ivar_pos
    //int ntrns;      // # of transport constituents on this material, in ivar_pos
    SIVAR_POSITION ivar_pos; // a list of indexes for independent variables on the element (not supermodel)
    int *ivar_loc;  // array of ivar_pos.n who's entries are location of independent variables (row id)
                    // element in **ivar array
    
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
void smat_physics_set_ivar_loc_array(SMAT_PHYSICS *m, int nmat_physics, SIVAR_POSITION* ivar_pos);
void smat_physics_set_ivar_pos(SMAT_PHYSICS *mat_physics_elem, int nmat_physics);
void smat_physics_set_nodal_pointers( SMAT_PHYSICS **mat_physics_node, SGRID *grid, int *elem1d_physics_mat,
    int *elem2d_physics_mat, int *elem3d_physics_mat, SMAT_PHYSICS *mat_physics_elem);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
