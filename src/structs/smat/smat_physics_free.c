/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_free.c This file frees an array of SMAT_PHYSICS structures          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat  (SMAT_PHYSICS *)  pointer to a transport material
 * @param[in]    nmat (int) number of materials to allocate
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_physics_free_array(SMAT_PHYSICS *mat, int nmat) {
    
    for (int imat=0; imat<nmat; imat++) {
        if(mat[imat].TRANSPORT!=NULL){
            mat[imat].TRANSPORT = (bool *) tl_free(sizeof(bool), mat[imat].ntrns, mat[imat].TRANSPORT);
        }
        if(mat[imat].model!=NULL){
            smodel_free_array(mat[imat].model, mat[imat].nSubmodels); 
        }
        if(mat[imat].ivar_loc != NULL){
            mat[imat].ivar_loc = (int *) tl_free(sizeof(bool), mat[imat].n, mat[imat].ivar_loc);
        }

    }
    mat = (SMAT_PHYSICS *) tl_free(sizeof(SMAT_PHYSICS), nmat, mat);
}