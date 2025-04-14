#ifndef _H_RESIDUAL_
#define _H_RESIDUAL_


void assemble_residual(SMODEL_SUPER *sm, SGRID *grid);
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int elem_nvars, int *elem_var_map, int imod, int nnodes, double scale);
void load_global_resid(double *residual, double *elem_rhs, int nnodes, int elem_nvars, int *local_dofs);

#endif
