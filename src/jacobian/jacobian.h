#ifndef _H_JACOBIAN_
#define _H_JACOBIAN_


void assemble_jacobian(SMODEL_SUPER *sm);
void load_elem_mat(SMODEL_SUPER *sm,  int elem_nnodes, int *nodeIDs, int mat_id, int elem_id, double **elem_mat);
void load_global_mat_split_CSR(double *vals, int *indptr, int *indices, double *off_diag_vals, int *off_diag_indptr, int *off_diag_indices, double **elem_mat, int ndofs_ele, int *dofs, int *global_dofs, int *local_range);
void load_global_mat_CSR(double *vals, int *indptr, int *indices, double **elem_mat, int ndofs_ele, int *global_dofs, int *local_range);
void perturb_var(double **elem_mat, SMODEL_SUPER *sm, SPDE *pde, 
    int ie, int mat_id, int nodes_on_element, int nvar_ele, int *elem_var_map, int npdes, int ele_var_no, int *NodeIDs, int DEBUG);
void elem_matrix_deriv(double **mat, int node_no, int var_no, int nnodes, int elem_nvars, double *local1, double *local2, double diff_ep);
#endif
