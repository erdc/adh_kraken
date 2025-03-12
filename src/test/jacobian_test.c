/*! \file jacobian_test.c This file tests the assembly of a Jacobian matrix */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests assembly of the Jacobian matrix
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int jacobian_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax) {
	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// return error code
    //++++++++++++++++++++++++++++++++++++++++++++++
	int ierr = -1;
	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Create a design model
    //++++++++++++++++++++++++++++++++++++++++++++++
	SMODEL_DESIGN dm; smodel_design_defaults(&dm);
	
	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Create a regular grid
    //++++++++++++++++++++++++++++++++++++++++++++++
	dm.grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
	//let's just do something simple
	//3x3 triangular element grid
	//double xmin = 0.0;
	//double xmax = 5.0;
	//int npx = 11;
	double theta = 0.0;
	double dz = 1.0;
	double a0 = -5.0;
	double ax = 0.0;
	double ax2 = 0.0;
	double ay = 0.0;
	double ay2 = 0.0;
	double axy = 0.0;
	double ax2y = 0.0;
	double axy2 = 0.0;
	double ax2y2 = 0.0;
	int flag3d =0;
    *(dm.grid) = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Reorder grid to minimize bandwidth
    //++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Fill in a simple design model
    //++++++++++++++++++++++++++++++++++++++++++++++
    double dt = 1.0;
	double t0 = 0.0;
	double tf = 1.0;
	int nSuperModels = 1;
	int *nphysics_mats;
	nphysics_mats = (int*) tl_alloc(sizeof(int), nSuperModels);
	nphysics_mats[0] = 1;
	//elemVarCode in general is triple pointer, nSuperModels x nphyics_mat[i] x 10
	char ***elemVarCode;
	elemVarCode = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	for (int i = 0; i< nSuperModels; i ++){
		elemVarCode[i] = (char **) tl_alloc(sizeof(char *), nphysics_mats[i]);
		for (int j = 0; j< nphysics_mats[i]; j++){
			elemVarCode[i][j] = (char *) tl_alloc(sizeof(char), 10);
		}
	}
	strcpy(&elemVarCode[0][0][0],"9"); // POISSON
	strcpy(&elemVarCode[0][0][1],"0"); // GW
	strcpy(&elemVarCode[0][0][2],"0"); // Transport


	//mat ids
	int **mat_ids;
	mat_ids = (int **) tl_alloc(sizeof(int*),nSuperModels);
	int nelems = dm.grid->nelems3d + dm.grid->nelems2d + dm.grid->nelems1d;
	for (int i = 0; i < nSuperModels; i++){
		mat_ids[i] = tl_alloc(sizeof(int), nelems);
		sarray_init_int(mat_ids[i],nelems);
	}

    smodel_design_init_no_read(&dm, dt, t0, tf, nSuperModels, nphysics_mats, elemVarCode, mat_ids);
    
    //sarray_init_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs);
    sarray_init_value_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs, 1.0);
    dm.superModel[0].LINEAR_PROBLEM = YES;


	assemble_jacobian(&(dm.superModel[0]));

	slin_sys_CSR_printScreen(dm.superModel[0].lin_sys);

	//free stuff
    smodel_design_free(&dm);
	return 0;
}
