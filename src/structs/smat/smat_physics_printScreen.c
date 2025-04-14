/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_printScreen.c This file prints an  SMAT_PHYSICS structures          */
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
 * @param[in] mat (SMAT_ELEM *)  pointer to a transport material
 * @param[in] id  (int) the SMAT ID
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_printScreen(SMAT_PHYSICS *m) {
    printf("----------------------------------------\n");
    printf("SMAT #%d || n: %d || ntrns: %d || npdes: %d\n",m->id+1,m->ivar_pos.n,m->ivar_pos.ntrns,m->npdes);
    printf("----------------------------------------\n");
    sivar_position_printScreen(&m->ivar_pos);
    for (int ipde=0; ipde<m->npdes; ipde++) {
        sadh_def_model_printScreen(&adh_def.model[m->pde[ipde].imod]);
    }
    printf("----------------------------------------\n");

}
