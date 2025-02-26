/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_update_array
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
void smat_physics_update_array(SMAT_PHYSICS *m, int nmat_physics, SIVAR_POSITION* ivar_pos){

	int ctr=0;
	int isubModel = 0;
	SMAT_PHYSICS *mat;
	//++++++++++++++++++++++++++++++++++++++++++++++
    // use the smodel_super level ivar_pos to
	// put correct indices in mat_physics
    //++++++++++++++++++++++++++++++++++++++++++++++
	for (int imat=0; imat<nmat_physics; imat++){
		//temporary variable, copy pointer
		mat = &(m[imat]);
		//check the IVAR_POSITION of this mat_physics
		if (mat->ivar_pos.h    != UNSET_INT) { mat->ivar_loc[ctr] = ivar_pos->h; ctr++;}
        if (mat->ivar_pos.u    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->u; ctr++;}
        if (mat->ivar_pos.v    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->v; ctr++;}
        if (mat->ivar_pos.w    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->w; ctr++;}
        if (mat->ivar_pos.uda  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->uda; ctr++;}
        if (mat->ivar_pos.vda  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->vda; ctr++;}
        if (mat->ivar_pos.dpl  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->dpl; ctr++;}
        if (mat->ivar_pos.prs  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->prs; ctr++;}     
        if (mat->ivar_pos.heat != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->heat; ctr++;}
        if (mat->ivar_pos.sal  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->sal; ctr++;}
        
        for (int nt = 0; nt<mat->ntrns ;nt++){
        	mat->ivar_loc[ctr] = ivar_pos->con[nt];
        	ctr++;
        }
        ctr=0;
        //++++++++++++++++++++++++++++++++++++++++++++++
    	// now use ivar_loc to inform the submodels
    	//++++++++++++++++++++++++++++++++++++++++++++++
    	// fill in model info
    	if (mat->SW1_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->h;
        	mat->model[isubModel].physics_vars[1] = ivar_pos->uda;
        	mat->model[isubModel].physics = SW1D_; // for body residuals
        	mat->model[isubModel].physics_init = SW1D_; 
        	isubModel++;
    	}
    	if (mat->SW2_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->h;
        	mat->model[isubModel].physics_vars[1] = ivar_pos->uda;
        	mat->model[isubModel].physics_vars[2] = ivar_pos->vda;
        	mat->model[isubModel].physics = SW2D_;
        	mat->model[isubModel].physics_init = SW2D_;
        	isubModel++;
    	}
    	if (mat->SW3_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->h;
        	mat->model[isubModel].physics_vars[1] = ivar_pos->u;
        	mat->model[isubModel].physics_vars[2] = ivar_pos->v;
        	mat->model[isubModel].physics = SW3D_;
        	mat->model[isubModel].physics_init = SW3D_;
        	isubModel++;
    	}
    	if (mat->NS3_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->u;
        	mat->model[isubModel].physics_vars[1] = ivar_pos->v;
        	mat->model[isubModel].physics_vars[2] = ivar_pos->w;
        	mat->model[isubModel].physics_vars[3] = ivar_pos->prs;
        	mat->model[isubModel].physics = NS3D_;
        	mat->model[isubModel].physics_init = NS3D_;
        	isubModel++;
    	}
    	if (mat->NS3_SPLIT) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->u;
        	mat->model[isubModel].physics_vars[1] = ivar_pos->v;
        	mat->model[isubModel].physics_vars[2] = ivar_pos->w;
        	mat->model[isubModel].physics = NS3DSPLIT_;
        	mat->model[isubModel].physics_init = NS3DSPLIT_;
        	isubModel++;
    	}
    	if (mat->DW_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->h;
        	mat->model[isubModel].physics = DW2D_;
        	mat->model[isubModel].physics_init = DW2D_;
        	isubModel++;
    	}
    	if (mat->WVEL_SPLIT) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->w;
        	mat->model[isubModel].physics = WVEL_;
        	mat->model[isubModel].physics_init = WVEL_;
        	isubModel++;
    	}
    	if (mat->PRESSURE) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->prs;
        	mat->model[isubModel].physics = PRS_;
        	mat->model[isubModel].physics_init = PRS_;
        	isubModel++;
    	}
    	if (mat->GW_FLOW) {
        	mat->model[isubModel].physics_vars[0] = ivar_pos->h;
        	mat->model[isubModel].physics = GW3D_;
        	mat->model[isubModel].physics_init = GW3D_;
        	isubModel++;
    	}
    	for (int itrns=0; itrns<mat->ivar_pos.ntrns; itrns++) {
        	if (mat->TRANSPORT[itrns]) {
            	mat->model[isubModel].physics_vars[0] = ivar_pos->con[itrns];
            	mat->model[isubModel].physics = TRNS_;
            	mat->model[isubModel].physics_init = TRNS_;
            	isubModel++;
        	}  
    	}
    	assert(isubModel == mat->nSubmodels);
    	isubModel=0;
	}

}