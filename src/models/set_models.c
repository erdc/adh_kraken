/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  set_models.c This file sets the global arrays of function pointers
 *  to all residual and init routines as created in include/model_codes.h             */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
void set_models(int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int), int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *)){
	//assign all values to residual vector of function pointers
	fe_resid[SW2D_] = fe_sw2_body_resid;
	fe_resid[POISSON2D_] = poisson_residual;
	fe_resid[HEAT2D_] = heat_residual;
	//assign all values to init vector of function pointers
	fe_init[SW2] = fe_sw2_init;

}
