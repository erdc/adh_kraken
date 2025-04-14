/*! \file  assemble_residual.c This file has functions responsible for assembling the global residual vector */
#include "adh.h"
static int DEBUG = OFF;
void load_pdes(SMODEL_SUPER *sm, int elem_nnodes, int *nodeIDs, int mat_id, int elem_id) {
    
    int i, imod = UNSET_INT;
    int var_code;
    int max_elem_dof = sm->ivar_pos.n*MAX_NNODE;
    int dofs[max_elem_dof];
    double elem_rhs[max_elem_dof];
    double eq_rhs[max_elem_dof];
    sarray_init_dbl(elem_rhs,max_elem_dof);
    SMAT_PHYSICS *mat = &(sm->mat_physics_elem[mat_id]); // alias
    int elem_nvars = mat->ivar_pos.n;
    int eq_vars[elem_nvars] ;
    for (int k=0;k<mat->npdes;k++){
        imod = mat->pde[k].imod;
        sarray_init_dbl(eq_rhs,max_elem_dof);
        //convention for filling temp will be:
        // for i in node (for j in nvar temp[nnode*i + j] = result)
        var_code = mat->pde[k].resid(sm,eq_rhs, elem_id, 0.0, UNSET_INT, UNSET_INT, 0, DEBUG);
        //add eq_rhs to elem_rhs
        //in order to do this we will need elemental vars and info about fe_resid routine
        //add_replace_elem_rhs(elem_rhs,eq_rhs,mat->pde[k].imod,elem_nnodes,-1.0);
        //Need to fix this, MARK, where are elem_nvars and elem_vars
        //what is relative position of equation variables??
        add_replace_elem_rhs(elem_rhs, eq_rhs, elem_nvars, mat->ivar_pos.var, imod, elem_nnodes, -1.0);
    }
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    //for residual we only need dof numbers local to process (including ghost nodes)
    //this is a complicated map but maybe we can simplify in simpler cases by replacing different routine
    //usually would take the local cell number and compute the associated dofs
    //but this has expanded arguments so it will work for elem1d,elem2d as well, cell # is implicit
    //++++++++++++++++++++++++++++++++++++++++++++++
    get_cell_dofs_ivars(dofs, sm->ivars, elem_nnodes, nodeIDs, mat->ivar_pos.n, mat->ivar_loc);

    //puts elem_rhs into global residual, applies Dirichlet conditions too?
    load_global_resid(sm->lin_sys->residual, elem_rhs, elem_nnodes, mat->ivar_pos.n, dofs);
}
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

    int max_elem_dof = sm->ivar_pos.n * MAX_NNODE;

    for (j=0;j<grid->nelems3d;j++){
        load_pdes(sm,grid->elem3d[j].nnodes,grid->elem3d[j].nodes,sm->elem3d_physics_mat[j],j);
    }
    for (j=0;j<grid->nelems2d;j++){
        load_pdes(sm,grid->elem2d[j].nnodes,grid->elem2d[j].nodes,sm->elem2d_physics_mat[j],j);
    }
    for (j=0;j<grid->nelems1d;j++){
        load_pdes(sm,grid->elem1d[j].nnodes,grid->elem1d[j].nodes,sm->elem1d_physics_mat[j],j);
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
void add_replace_elem_rhs(double *elem_rhs, double *eq_rhs, int elem_nvars, int *elem_var_map, int imod, int nnodes, double scale){
    //for each node, place the rhs entries of a specific pde residual routine in the correct slots of an elemental rhs

    int inode,eq_var; //eq_nvars;
    int current_var;
    //number of digits in eq_vars will be the number of variables in this residual
    //eq_nvars = count_digits(eq_vars);
    int neq_vars = adh_def.model[imod].nivars;

    for (inode=0;inode<nnodes;inode++){

        for (eq_var=0;eq_var<neq_vars;eq_var++){

            //get var index within elem vars
            //current_var = eq_vars[eq_var];
            //get var index within elem vars
            current_var = elem_var_map[adh_def.model[imod].var[eq_var]];

            //maps the current var from the residual to the correct var number in elem_vars

            //now we know what the current_var is in the whole element, put those entries into the elem_rhs
            elem_rhs[inode*elem_nvars + current_var] += scale*eq_rhs[inode*neq_vars+eq_var];

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
