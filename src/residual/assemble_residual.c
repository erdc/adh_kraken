/*! \file  assemble_residual.c This file has functions responsible for assembling the global residual vector */
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function assembles the global residual vector using elemental resid routines
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] sm (SMODEL_SUPER*) - the super model where the residual resides
 *  @param[in] grid (SGRID*) - the grid over which the monolithic residual resides
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void assemble_residual(SMODEL_SUPER *sm, SGRID *grid) {

    int j,k;
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // zero out residual
    // create small temporary arrays
    //++++++++++++++++++++++++++++++++++++++++++++++

    sarray_init_dbl(sm->lin_sys->residual, (sm->ndofs));

    double elem_rhs[MAX_ELEM_DOF];
    double eq_rhs[MAX_ELEM_DOF];
    int dofs[MAX_ELEM_DOF];
    int nnodes;
    int physics_vars[MAX_NVAR];
    int var_code;
    int nvars_elem, nphysics_models, mat_id, nvar_pde;
    int **ivars = sm->ivars;

    int offset;

    //loop through and assemble residual
    for (j=0;j<grid->nelems3d;j++){
        sarray_init_dbl(elem_rhs,MAX_ELEM_DOF);

        nnodes = grid->elem3d[j].nnodes;
        mat_id = sm->elem3d_physics_mat[j];
        nvars_elem =  sm->mat_physics_elem[mat_id].n;
        nphysics_models = sm->mat_physics_elem[mat_id].nSubmodels;
        offset = sm->resid_ptr[mat_id];


        for (k=0;k<nphysics_models;k++){
            sarray_init_dbl(eq_rhs,MAX_ELEM_DOF);
            
            nvar_pde = sm->mat_physics_elem[mat_id].model[k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->mat_physics_elem[mat_id].model[k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)
            //var_code = smodel_super_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG, fe_resid[sm->mat_physics_elem[mat_id].model[k].physics]);
            //replace with:
            var_code = sm->fe_resid[mat_id+k](sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            //var_code = smodel_super_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG, fe_resid[sm->mat_physics_elem[mat_id].model[k].physics]);
            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,nvar_pde,physics_vars,nnodes,-1.0);
        }
        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem3d[j].nodes, nvars_elem, sm->mat_physics_elem[mat_id].ivar_loc);

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->lin_sys->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }

    for (j=0;j<grid->nelems2d;j++){

        sarray_init_dbl(elem_rhs,MAX_ELEM_DOF);
        nnodes = grid->elem2d[j].nnodes;
        mat_id = sm->elem2d_physics_mat[j];
        nvars_elem =  sm->mat_physics_elem[mat_id].n;
        nphysics_models = sm->mat_physics_elem[mat_id].nSubmodels;


        for (k=0;k<nphysics_models;k++){

            sarray_init_dbl(eq_rhs,MAX_ELEM_DOF);
            
            nvar_pde = sm->mat_physics_elem[mat_id].model[k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->mat_physics_elem[mat_id].model[k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)
            //var_code = smodel_super_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG, fe_resid[sm->mat_physics_elem[mat_id].model[k].physics]);
            printf("In assemble_residual\n");
            var_code = sm->fe_resid[mat_id+k](sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,nvar_pde,physics_vars,nnodes,-1.0);
        }
        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem3d[j].nodes, nvars_elem, sm->mat_physics_elem[mat_id].ivar_loc);

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->lin_sys->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }

    for (j=0;j<grid->nelems1d;j++){

        sarray_init_dbl(elem_rhs,MAX_ELEM_DOF);
        nnodes = grid->elem1d[j].nnodes;
        mat_id = sm->elem1d_physics_mat[j];
        nvars_elem =  sm->mat_physics_elem[mat_id].n;
        nphysics_models = sm->mat_physics_elem[mat_id].nSubmodels;


        for (k=0;k<nphysics_models;k++){
            sarray_init_dbl(eq_rhs,MAX_ELEM_DOF);
            
            nvar_pde = sm->mat_physics_elem[mat_id].model[k].nvar;
            sarray_init_int(physics_vars, nvar_pde);
            sarray_copy_int(physics_vars, sm->mat_physics_elem[mat_id].model[k].physics_vars,nvar_pde);
            //convention for filling temp will be:
            // for i in node (for j in nvar temp[nnode*i + j] = result)
            //var_code = smodel_super_resid(sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG, fe_resid[sm->mat_physics_elem[mat_id].model[k].physics]);
            var_code = sm->fe_resid[mat_id+k](sm,eq_rhs,j, 0.0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
            //add eq_rhs to elem_rhs
            //in order to do this we will need elemental vars and info about fe_resid routine
            add_replace_elem_rhs(elem_rhs,eq_rhs,nvars_elem,nvar_pde,physics_vars,nnodes,-1.0);
        }
        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        //for residual we only need dof numbers local to process (including ghost nodes)
        //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
        //usually would take the local cell number and compute the associated dofs
        //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem3d[j].nodes, nvars_elem, sm->mat_physics_elem[mat_id].ivar_loc);

        //puts elem_rhs into global residual, applies Dirichlet conditions too?
        load_global_resid(sm->lin_sys->residual, elem_rhs, nnodes, nvars_elem, dofs);
    }

   
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function adds the residual of one PDE into the elemental residual
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] elem_rhs (double*) - array containing the elemental residual
 *  @param[in] eq_rhs (double*) - array containg the residual from a single PDE routine
 *  @param[in] elem_nvars (int) - integer that is the length of *elem_nvars
 *  @param[in] elem_vars (int*) - array containg the unique variable codes defined on an element
 *  @param[in] eq_nvars (int) - integer that is the length of *eq_vars
 *  @param[in] eq_vars (int*) - array containg the unique variable codes defined by the PDE
 *  @param[in] nnodes (int) - integer that is the number of nodes on one element
 *  @param[in] scale (double) - scalar factor to scale the residual
 * 
 *  \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int elem_nvars, int eq_nvars, int *eq_vars, int nnodes, double scale){

    //for each node, place the rhs entries of a specific pde residual routine in the correct slots of an elemental rhs

    int inode,eq_var; //eq_nvars;
    int current_var;
    bool notFound = TRUE;
    int k,save_k;
    //number of digits in eq_vars will be the number of variables in this residual
    //eq_nvars = count_digits(eq_vars);

    for (inode=0;inode<nnodes;inode++){

        for (eq_var=0;eq_var<eq_nvars;eq_var++){

            //get var index within elem vars
            current_var = eq_vars[eq_var];

            //maps the current var from the residual to the correct var number in elem_vars

            //now we know what the current_var is in the whole element, put those entries into the elem_rhs
            elem_rhs[inode*elem_nvars + current_var] += scale*eq_rhs[inode*eq_nvars+eq_var];

        }
    
    } 

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes an elemental residual and loads into the global residual vector
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] residual (double*) - the global residual
 *  @param[in] elem_rhs (double*) - the local element right-hand-side
 *  @param[in] nnodes (int) - the number of local nodes on the element
 *  @param[in] elem_nvars (int) - the number of active variables on the element
 *  @param[in] local_dofs (int*) - array of the degree of freedom numbers local to the processor
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void load_global_resid(double *residual, double *elem_rhs, int nnodes, int elem_nvars, int *local_dofs) {
    int index,i;
    /// assembles global residual
    for (i=0; i<nnodes*elem_nvars; i++) {

            //map the current var from the residual to the correct var number in processor residual
            //using local_dofs map

            index = local_dofs[i];
            //does minus convention hold for all residuals?
            residual[index]     += elem_rhs[i];
                    
        }

}
