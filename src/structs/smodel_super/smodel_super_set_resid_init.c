/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super_set_resid_init
 * The file is for filling in the ivar_positions in an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and sets up physics material properties
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
void smodel_super_set_resid_init(SMODEL_SUPER *sm, SMAT_PHYSICS *mat_physics, int nmat_physics){
	int n_sub_mods;
	int ctr = 0;

	sm->resid_ptr = (int *) tl_alloc(sizeof(int),nmat_physics);
	sm->init_ptr = (int *) tl_alloc(sizeof(int),nmat_physics);
	sarray_init_int(sm->resid_ptr, nmat_physics); 
	sarray_init_int(sm->init_ptr, nmat_physics); 
	// Loop over each mat physics and store the proper function pointers
	// canonical ordering is flattened version of
	// loop over mat physics then each submodel in that physics
	for (int imat=0; imat<nmat_physics; imat++){
		n_sub_mods= mat_physics[imat].nSubmodels;
		sm->resid_ptr[imat] = ctr;
		//redundant? what about resids that share init routine??
		sm->init_ptr[imat] = ctr;
		for(int nmod=0; nmod<n_sub_mods ; nmod++){
			sm->fe_resid[ctr] = &mat_physics[imat].model[nmod].physics;
			sm->fe_init[ctr] = &mat_physics[imat].model[nmod].physics_init;
			ctr++;
		}
	}


}
