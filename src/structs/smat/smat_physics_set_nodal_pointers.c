/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_set_nodal_pointers.c The file is for allocating an initializing an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 * @param[in]    nmat        (int) number of materials to allocate
 * @param[in]    codes       (char [][10]) an array of strings that encode the materials physics
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_set_nodal_pointers( SMAT_PHYSICS **mat_physics_node, SGRID *grid, int *elem1d_physics_mat,
    int *elem2d_physics_mat, int *elem3d_physics_mat, SMAT_PHYSICS *mat_physics_elem){

    int ie,imat,nd,i;
    int *node_nvars = tl_alloc(sizeof(int), grid->nnodes); 
    sarray_init_int(node_nvars,grid->nnodes);

    for (ie=0; ie<grid->nelems1d; ie++) {
        if (grid->elem1d[ie].bflag != BODY) continue;
        imat = elem1d_physics_mat[ie];
        for (i=0; i<grid->elem1d[ie].nnodes; i++) {
            nd = grid->elem1d[ie].nodes[i];
            if (node_nvars[nd] < mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = mat_physics_elem[imat].ivar_pos.n;
                mat_physics_node[nd] =  &(mat_physics_elem[elem1d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems2d; ie++) {
        if (grid->elem2d[ie].bflag != BODY) continue;
        imat = elem2d_physics_mat[ie];
        for (i=0; i<grid->elem2d[ie].nnodes; i++) {
            nd = grid->elem2d[ie].nodes[i];
            if (node_nvars[nd] < mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = mat_physics_elem[imat].ivar_pos.n;
                mat_physics_node[nd] =  &(mat_physics_elem[elem2d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems3d; ie++) {
        imat = elem3d_physics_mat[ie];
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            nd = grid->elem3d[ie].nodes[i];
            if (node_nvars[nd] < mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = mat_physics_elem[imat].ivar_pos.n;
                mat_physics_node[nd] =  &(mat_physics_elem[elem3d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }

    node_nvars = tl_free(sizeof(int), grid->nnodes,node_nvars); 
}