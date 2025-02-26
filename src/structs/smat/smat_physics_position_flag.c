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

    sarray_init_int(FLAG,MAX_TRNS_VARS + MAX_VARS);
    
    for (int i=0; i<nnodes; i++) {
        //smat_elem_printScreen(mat_node[i],i);
        if (mat_node[i]->ivar_pos.h    != UNSET_INT) FLAG[0] = 1;
        if (mat_node[i]->ivar_pos.u    != UNSET_INT) FLAG[1] = 1;
        if (mat_node[i]->ivar_pos.v    != UNSET_INT) FLAG[2] = 1;
        if (mat_node[i]->ivar_pos.w    != UNSET_INT) FLAG[3] = 1;
        if (mat_node[i]->ivar_pos.uda  != UNSET_INT) FLAG[4] = 1;
        if (mat_node[i]->ivar_pos.vda  != UNSET_INT) FLAG[5] = 1;
        if (mat_node[i]->ivar_pos.dpl  != UNSET_INT) FLAG[6] = 1;
        if (mat_node[i]->ivar_pos.prs  != UNSET_INT) FLAG[7] = 1;       
        if (mat_node[i]->ivar_pos.heat != UNSET_INT) FLAG[8] = 1;
        if (mat_node[i]->ivar_pos.sal  != UNSET_INT) FLAG[9] = 1;

        for (int itrns=0; itrns<MAX_TRNS_VARS ; itrns++) {
            // CJT - note: Different types of transport may cause issues here
            if (mat_node[i]->ivar_pos.con[itrns] != UNSET_INT) {
                FLAG[MAX_VARS + itrns] = 1;   
            }
        }

    }
}