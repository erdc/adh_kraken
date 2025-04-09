/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  set_function_pointers.c This file sets the global arrays of function pointers
 *  to all residual and init routines as created in include/model_codes.h             */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"


void set_function_pointers(  int (*forward_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*)){
	

    // we can actually add empty routine as placeholder
    // add actual routines as they come
//    resid_routines[NO_RESID_] = no_resid;
//    resid_routines[SW1D_] = no_resid;
//	resid_routines[SW2D_] = fe_sw2_body_resid;
//	resid_routines[SW3D_] = no_resid;
//	resid_routines[WVEL_] = no_resid;
//	resid_routines[NS1D_] = no_resid;
//	resid_routines[NS2D_] = no_resid;
//	resid_routines[NS3D_] = no_resid;
//	resid_routines[NS1DSPLIT_] = no_resid;
//	resid_routines[NS2DSPLIT_] = no_resid;
//	resid_routines[NS3DSPLIT_] = no_resid;
//	resid_routines[PRS_] = no_resid;
//	resid_routines[DW1D_] = no_resid;
//	resid_routines[DW2D_] = no_resid;
//	resid_routines[GW1D_] = no_resid;
//	resid_routines[GW2D_] = no_resid;
//	resid_routines[GW3D_] = no_resid;
//	resid_routines[TRNS_] = no_resid;
//	resid_routines[TRNS1D_] = no_resid;
//	resid_routines[TRNS2D_] = no_resid;
//	resid_routines[TRNS3D_] = no_resid;
//	resid_routines[HEAT1D_] = no_resid;
//	resid_routines[HEAT2D_] = heat_residual;
//	resid_routines[HEAT3D_] = no_resid;
//	resid_routines[POISSON2D_] = poisson_residual;
//	// Definitions for weak boundary treatments, included as resid routines
//	resid_routines[SW2_BC_DISCHARGE_] =  fe_sw2_bc_discharge;
//	resid_routines[SW2_BC_ELE_] = fe_sw2_bc_ele;
//	resid_routines[SW2_BC_FLAPD_] = fe_sw2_bc_flapd;
//	resid_routines[SW2_BC_FLAPU_] = fe_sw2_bc_flapu;
//	resid_routines[SW2_BC_FLUX_] = fe_sw2_bc_flux;
//	resid_routines[SW2_BC_H_] = fe_sw2_bc_h;
//	resid_routines[SW2_BC_HYBRID_] = fe_sw2_bc_hybrid;
//	resid_routines[SW2_BC_OUTFLOW_] = fe_sw2_bc_outflow;
//	resid_routines[SW2_BC_SLUICED_] = fe_sw2_bc_sluiced;
//	resid_routines[SW2_BC_SLUICEU_] = fe_sw2_bc_sluiceu;
//	resid_routines[SW2_BC_TAILWATER_] = fe_sw2_bc_tailwater;
//	resid_routines[SW2_BC_VEL_] = fe_sw2_bc_vel;
//	resid_routines[SW2_BC_VEL_ELE_] = fe_sw2_bc_vel_ele;
//	resid_routines[SW2_BC_WEIRD_] = fe_sw2_bc_weird;
//	resid_routines[SW2_BC_WEIRU_] = fe_sw2_bc_weiru;


	//assign all values to init vector of function pointers
//	init_routines[NO_INIT] = no_init;
//	init_routines[SW2_INIT] = fe_sw2_init;

	//assign all pointers to time loop
	//so far only one but maybe RK in future
	forward_stepper[FE_NEWTON] = fe_newton;

}
