#ifndef H_SDVAR_
#define H_SDVAR_

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/

typedef struct {
    int n;            // the number of active dependent variables
    int  *var;        // position in array form
    bool *print_flag; // true if the variable is output
} SDVAR_POSITION;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
typedef struct {

    // node-based variables
    int n_dvar;     // number of active dependent variables
    int nnode_dvar; // number of active nodes
    SDVAR_POSITION sdvar_pos_node; // variable positions in matrix
    double **nodal_dvar; // [n_dvar][nnode_dvar]

    // element-based surface water variables
    int n_dvar_elem_dbl;
    int n_dvar_elem_int;
    int nelem_dvar;//number of active elements
    SDVAR_POSITION sdvar_pos_elem_dbl;
    SDVAR_POSITION sdvar_pos_elem_int;
    double **elem_dvar;//[n_dvar_elem_dbl][nelem_dvar]
    int **elem_flags;//[n_dvar_elem_int][nelem_dvar]

    //maps
    int *dvar_node_map; //[nnode] array goes NodeID->index within nodal_dvar
    int *dvar_active; //n_dvar array that has the NodeID, convenient for printing out stuff
    int *dvar_elem_map; //nelem array takes elem # -> index within elem_var

} SDVAR;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/

void sdvar_position_init(SDVAR_POSITION *dp);
void sdvar_position_free(SDVAR_POSITION *dp);
void sdvar_position_printScreen(SDVAR_POSITION *dp);
void sdvar_alloc_init(SDVAR *sdvar, int nnode, int n_dvar, int n_dvar_elem_dbl, int n_dvar_elem_int, int nnode_dvar, int nelem_dvar);
#endif

