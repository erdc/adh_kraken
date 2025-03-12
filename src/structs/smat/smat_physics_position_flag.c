/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_position_flag.c This file flags variable positions on the element   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Flags independent variables on all nodes to get all ivariables used on superModel
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in]    mat (SMAT_ELEM *)  pointer to a transport material
 * @param[in]    nnodes (int) the number of nodes on the grid
 * @param[inout] FLAG (int *) flags which variables are being used over the whole set of nodes
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_position_flag(SMAT_PHYSICS **mat_node, int nnodes, int *FLAG) {
    sarray_init_int(FLAG,N_IVARS_TOTAL);
    for (int i=0; i<nnodes; i++) {
        for (int ivar=0; ivar<N_IVARS_TOTAL; ivar++) {
            //printf("mat_node[%d]->ivars.var[ivar=%d]: %d\n",i,ivar,mat_node[i]->ivars.var[ivar]);
            if (mat_node[i]->ivar_pos.var[ivar] != UNSET_INT) FLAG[ivar] = 1;
        }
    }
}
