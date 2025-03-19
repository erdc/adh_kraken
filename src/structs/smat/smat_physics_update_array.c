/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_update_array
 * The file is for filling in the ivar_positions in an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Sets up the ivar_loc array within smat_physics, necessary to map elemental
 *  matrices to global
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
void smat_physics_update_array(SMAT_PHYSICS *m, int nmat_physics, SIVAR_POSITION* ivar_pos){
	int ctr=0;
	SMAT_PHYSICS *mat;
	//++++++++++++++++++++++++++++++++++++++++++++++
    // use the smodel_super level ivar_pos to
	// put correct indices in mat_physics
    //++++++++++++++++++++++++++++++++++++++++++++++
	for (int imat=0; imat<nmat_physics; imat++){

		//temporary variable, copy pointer
		mat = &(m[imat]);
        if (mat->ivar_pos.var[_H]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_H]; ctr++;}
        if (mat->ivar_pos.var[_U]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_U]; ctr++;}
        if (mat->ivar_pos.var[_V]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_V]; ctr++;}
        if (mat->ivar_pos.var[_W]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_W]; ctr++;}
        if (mat->ivar_pos.var[_UDA]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_UDA]; ctr++;}
        if (mat->ivar_pos.var[_VDA]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_VDA]; ctr++;}
        if (mat->ivar_pos.var[_DPL]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_DPL]; ctr++;}
        if (mat->ivar_pos.var[_PRS]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_PRS]; ctr++;}     
        if (mat->ivar_pos.var[_HEAT] != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_HEAT]; ctr++;}
        if (mat->ivar_pos.var[_SAL]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_SAL]; ctr++;}   
        for (int nt = 0; nt<mat->ntrns ;nt++){
            mat->ivar_loc[ctr] = ivar_pos->var[N_IVARS + nt];
            ctr++;
        }

        ctr=0;
        
	}

}