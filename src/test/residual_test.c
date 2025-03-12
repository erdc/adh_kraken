/*! \file residual_test.c This file tests the residual assembly */
#include "adh.h"
static void find_analytic_residual_linear_poisson(double *resid, SGRID *grid, double hx, double hy);
static double RESID_TEST_TOL = 1e-12;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the assembly of a residual vector
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int residual_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax) {
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

    sgrid_reorder(dm.grid,2);

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
    
    sarray_init_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs);
    dm.superModel[0].LINEAR_PROBLEM = YES;


	assemble_residual(&(dm.superModel[0]), dm.grid);

	//print final residual
	//sarray_printScreen_dbl(dm.superModel[0].lin_sys->residual, dm.superModel[0].ndofs, "residual");

	double *exact_sol;
	int nnodes = dm.grid->nnodes;
	double hx = (xmax - xmin)/(npx-1);
	double hy = (ymax - ymin)/(npy-1);
	exact_sol = (double *) tl_alloc(sizeof(double), nnodes);
	sarray_init_dbl(exact_sol, nnodes);

	find_analytic_residual_linear_poisson(exact_sol, dm.grid, hx, hy);

	//compute L2 and Linf error
	double l2_err =  l2_error(dm.superModel[0].lin_sys->residual, exact_sol, nnodes);
	double linf_err =  linf_error(dm.superModel[0].lin_sys->residual, exact_sol, nnodes);



	//return -1 if failed, 0 if good
	
	if(l2_err < RESID_TEST_TOL && linf_err < RESID_TEST_TOL){
		ierr=0;
	}


    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Free Memory
    //++++++++++++++++++++++++++++++++++++++++++++++
    exact_sol = tl_free(sizeof(double), nnodes, exact_sol);
    smodel_design_free(&dm);
    
#ifdef _DEBUG
    printf(">assemble residual error code %d : l2err = %16.4e, linf err=%16.4e \n", ierr, l2_err, linf_err);
#endif
	return ierr;
}



void find_analytic_residual_linear_poisson(double *resid, SGRID *grid, double hx, double hy){
	//computes the integral int_f_v_dx where f=-6
	//h is the element width

	//this will be (n elements connected to node) * -6 * (h^2)/6

	//get connections per node
	int nd1,nd2,nd3;
	int *n_connections;
	int nnodes = grid->nnodes;
	n_connections = (int *) tl_alloc(sizeof(int), nnodes);
	sarray_init_int(n_connections, nnodes);

	for (int ie = 0; ie<grid->nelems2d ; ie++ ){
		nd1 = grid->elem2d[ie].nodes[0];
		nd2 = grid->elem2d[ie].nodes[1];
		nd3 = grid->elem2d[ie].nodes[2];

		n_connections[nd1]++;
		n_connections[nd2]++;
		n_connections[nd3]++;
	} 


	for (int i = 0; i<nnodes; i++){
		resid[i] = n_connections[i]*(-6.0)*hx*hy/6.0;
	}

	return;

}