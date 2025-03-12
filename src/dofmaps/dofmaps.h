#ifndef H_DOFMAPS_
#define H_DOFMAPS_


//cg_maps

//uses array look up

void get_cell_dofs_ivars(int *local_dofs, int **ivars, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_var_pos);
void global_to_local_dbl_ivars(double *local_vals, int *nodeIDs, int nnodes, int *map_array, double *global_vals);
void global_to_local_SVECT2D_ivars(SVECT2D *local_vals, int *nodeIDs, int nnodes, int **ivars, int varx, int vary, double *global_vals);

//general_maps
void local_dofs_to_global_dofs(int *global_dofs,int ndofs_on_ele,int *dofs,int *local_range,int *ghosts);
int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost);


#endif
