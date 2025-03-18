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
    int h_loc, u_loc, v_loc, w_loc, uda_loc, vda_loc, dpl_loc, prs_loc, heat_loc, sal_loc;
    int con_loc[N_IVARS_TRNS];
	//++++++++++++++++++++++++++++++++++++++++++++++
    // use the smodel_super level ivar_pos to
	// put correct indices in mat_physics
    //++++++++++++++++++++++++++++++++++++++++++++++
	for (int imat=0; imat<nmat_physics; imat++){

		//temporary variable, copy pointer
		mat = &(m[imat]);
		//check the IVAR_POSITION of this mat_physics
		if (mat->ivar_pos.var[_H]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_H]; h_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_U]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_U]; u_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_V]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_V]; v_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_W]    != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_W]; w_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_UDA]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_UDA]; uda_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_VDA]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_VDA]; vda_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_DPL]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_DPL]; dpl_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_PRS]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_PRS]; prs_loc=ctr; ctr++;}     
        if (mat->ivar_pos.var[_HEAT] != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_HEAT]; heat_loc=ctr; ctr++;}
        if (mat->ivar_pos.var[_SAL]  != UNSET_INT) {mat->ivar_loc[ctr] = ivar_pos->var[_SAL]; sal_loc=ctr; ctr++;}
        
        for (int nt = 0; nt<mat->ntrns ;nt++){
        	mat->ivar_loc[ctr] = ivar_pos->var[N_IVARS + nt];
            con_loc[nt] = ctr;
        	ctr++;
        }

        ctr=0;
        //++++++++++++++++++++++++++++++++++++++++++++++
    	// now use ivar_loc to inform the submodels
    	//++++++++++++++++++++++++++++++++++++++++++++++
    	// fill in model info
    	if (mat->SW1_FLOW) {
        	mat->model[isubModel].physics_vars[0] = h_loc;
        	mat->model[isubModel].physics_vars[1] = uda_loc;
        	mat->model[isubModel].physics = SW1D_; // for body residuals
        	mat->model[isubModel].physics_init = SW1D_; 
        	isubModel++;
    	}
    	if (mat->SW2_FLOW) {
        	mat->model[isubModel].physics_vars[0] = h_loc;
        	mat->model[isubModel].physics_vars[1] = uda_loc;
        	mat->model[isubModel].physics_vars[2] = vda_loc;
        	mat->model[isubModel].physics = SW2D_;
        	mat->model[isubModel].physics_init = SW2_INIT;
        	isubModel++;
    	}
    	if (mat->SW3_FLOW) {
        	mat->model[isubModel].physics_vars[0] = h_loc;
        	mat->model[isubModel].physics_vars[1] = u_loc;
        	mat->model[isubModel].physics_vars[2] = v_loc;
        	mat->model[isubModel].physics = SW3D_;
        	mat->model[isubModel].physics_init = SW3D_;
        	isubModel++;
    	}

    	if (mat->NS3_FLOW) {
        	mat->model[isubModel].physics_vars[0] = u_loc;
        	mat->model[isubModel].physics_vars[1] = v_loc;
        	mat->model[isubModel].physics_vars[2] = w_loc;
        	mat->model[isubModel].physics_vars[3] = prs_loc;
        	mat->model[isubModel].physics = NS3D_;
        	mat->model[isubModel].physics_init = NS3D_;
        	isubModel++;
    	}
    	if (mat->NS3_SPLIT) {
        	mat->model[isubModel].physics_vars[0] = u_loc;
        	mat->model[isubModel].physics_vars[1] = v_loc;
        	mat->model[isubModel].physics_vars[2] = w_loc;
        	mat->model[isubModel].physics = NS3DSPLIT_;
        	mat->model[isubModel].physics_init = NS3DSPLIT_;
        	isubModel++;
    	}
    	if (mat->DW_FLOW) {

        	mat->model[isubModel].physics_vars[0] = h_loc;
        	mat->model[isubModel].physics = DW2D_;
        	mat->model[isubModel].physics_init = DW2D_;
        	isubModel++;
    	}
    	if (mat->WVEL_SPLIT) {
        	mat->model[isubModel].physics_vars[0] = w_loc;
        	mat->model[isubModel].physics = WVEL_;
        	mat->model[isubModel].physics_init = WVEL_;
        	isubModel++;
    	}
    	if (mat->PRESSURE) {
        	mat->model[isubModel].physics_vars[0] = prs_loc;
        	mat->model[isubModel].physics = PRS_;
        	mat->model[isubModel].physics_init = PRS_;
        	isubModel++;
    	}
    	if (mat->GW_FLOW) {
        	mat->model[isubModel].physics_vars[0] = h_loc;
        	mat->model[isubModel].physics = GW3D_;
        	mat->model[isubModel].physics_init = GW3D_;
        	isubModel++;
    	}
        if (mat->POISSON) {
            mat->model[isubModel].physics_vars[0] = h_loc;
            mat->model[isubModel].physics = POISSON2D_;
            //leave unset
            //mat->model[isubModel].physics_init = NO_INIT_;
            isubModel++;
        }
        if (mat->HEAT) {
            mat->model[isubModel].physics_vars[0] = h_loc;
            mat->model[isubModel].physics = HEAT2D_;
            //leave unset
            //mat->model[isubModel].physics_init = NO_INIT_;
            isubModel++;
        }

    	for (int itrns=0; itrns<mat->ivar_pos.ntrns; itrns++) {
        	if (mat->TRANSPORT[itrns]) {
            	mat->model[isubModel].physics_vars[0] = con_loc[itrns];
            	mat->model[isubModel].physics = TRNS_;
            	mat->model[isubModel].physics_init = TRNS_;
            	isubModel++;
        	}  
    	}


    	assert(isubModel == mat->nSubmodels);
    	isubModel=0;
	}

}