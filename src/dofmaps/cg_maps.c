/*! \file  cg_maps.c This file collections functions responsible for finding order of finite element
 * degrees of freedom for CG elements  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Routine that gives an array of degrees of freedom local to the current process for a CG element using ivars double array
 *             Designed to work for nodal (CG) based dof mappings
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_dofs (int*) - an array of integers that will give the degree of freedom numbers (equation numbers) for a
 *  given element local to the process
 *  @param[in] fmaplocal (int*) - an array of integers that gives the lowest d.o.f at a given node
 *  @param[in] nnodes (int) - the number of nodes on the element
 *  @param[in] local_node_ids (int*) - array of length nnodes containing node numbers local to process
 *  @param[in] elem_nvars (int) - number of solution variables active on the element
 *  @param[in] elem_vars (int*) - array of length elem_nvars that has the integer code for each variable
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - an array of SMAT_PHYSICS structs that contains variable info for each node
 *  @param[in] nodal_physics_mat_id (int*) - array of integers that gives the nodal physics mat id
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void get_cell_dofs_ivars(int *local_dofs, int **ivars, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_var_pos){
    int i,j,ctr, nodeID;
    int current_var;
    ctr = 0;
    //ivars will only work for CG, need to rethink for DG or possibly mixed CG-DG materials
    for (i=0; i<nnodes; i++){
        nodeID=local_node_ids[i];
        for (j=0; j<elem_nvars; j++){
            //map the current var from the residual to the correct var number in global residual
            current_var = elem_var_pos[j];
            //loop through the nodal vars to look for match
#ifdef _DEBUG
            assert(ivars[current_var][nodeID] != UNSET_INT);
#endif
            local_dofs[ctr] =  ivars[current_var][nodeID];
            ctr++;
         }
     }

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function takes a global array of values, and picks out the values for an array of nodes
 *  given a map array which takes nodeID -> dof index for the specific variable. (works for CG only)
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] local_vals (double*) - array of values corresponding to the node Ids
 *  @param[in] nodeIDs (int*) - node IDs (from the grid)
 *  @param[in] nnodes (int) - number of nodes requested to pull out from global array
 *  @param[in] map_array (int*) -  array that takes node ID and produces dof # for a specific variable
 *  @param[in] global_vals (double*) -  the global array of values that we want to pull from
 * \note 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_dbl_ivars(double *local_vals, int *nodeIDs, int nnodes, int *map_array, double *global_vals){
        for(int i =0;i<nnodes;i++){
                local_vals[i] = global_vals[map_array[nodeIDs[i]]];
        }
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     extracts sub-array of solution values for specific variable at a given set of nodes
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local (double*) - the array containing local (to the element) values of a specific variable
 *  @param[in] global (double*) - the full array of the solution vector
 *  @param[in] nodes (int*) - array of node IDs local to process
 *  @param[in] nnodes (int) - the numer of nodes in the nodes array
 *  @param[in] var (int) - the variable code to be extracted
 *  @param[in] fmaplocal (int*) - array containing first dof number (local to process) at each node
 *  @param[in] node_physics_mat (SMAT_PHYSICS*) - array of SMAT_PHYSICS structs containing variable info at nodes
 *  @param[in] nodal_physics_mat_id (int*) - array of ints containing the nodal physics mat id at each node
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void global_to_local_SVECT2D_ivars(SVECT2D *local_vals, int *nodeIDs, int nnodes, int **ivars, int varx, int vary, double *global_vals) {
    int i=0;
    int temp1,temp2;
    for (i=0; i<nnodes; i++) {
        local_vals[i].x = global_vals[ ivars[varx][nodeIDs[i]]   ];
        local_vals[i].y = global_vals[ ivars[vary][nodeIDs[i]]   ];
    }
}




