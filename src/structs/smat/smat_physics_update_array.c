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
        for (int var= 0 ; var<N_IVARS_TOTAL; var++){
            //if variable is active on element, get its position in **ivars
            if (mat->ivar_pos.var[var]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[var]; ctr++;}
        }
        ctr=0;
	}

}