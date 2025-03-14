/*! \file newton_test.c This file tests the PETSc solver for split CSR matrix */
#include "adh.h"
static double LINEAR_NEWTON_TEST_TOL = 1e-7;
static double NONLINEAR_NEWTON_TEST_TOL = 1e-3;
static int linear_newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
static void compute_exact_solution_poisson(double *u_exact, int ndof, SGRID *grid);
static int nonlinear_newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
static void compute_exact_solution_nonlinear_poisson(double *u_exact, int ndof, SGRID *grid);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the Newton solvet using a Poisson equation with analytic
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax){
	int err = linear_newton_test(npx, npy, xmin, xmax, ymin, ymax);
	err+= nonlinear_newton_test(npx, npy, xmin, xmax, ymin, ymax);
	return err;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the Newton solvet using a Poisson equation with analytic
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int linear_newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax) {
	
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
	int err_code=-1;
    *(dm.grid) = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );
    int nnodes;
    nnodes = dm.grid->nnodes;
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
    
    //sarray_init_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs);
    sarray_init_value_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs, 1.0);
    dm.superModel[0].LINEAR_PROBLEM = YES;

 	SMODEL_SUPER *sm;
	sm = &(dm.superModel[0]);

	printf("SETTING UP BCMASK\n");

	// intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<sm->ndofs; local_index++){

		dm.superModel[0].dirichlet_data[local_index] = 0.0;
		dm.superModel[0].sol_old[local_index] = 20.0;
		dm.superModel[0].sol[local_index] = 20.0;
		dm.superModel[0].lin_sys->dsol[local_index] = 0.0;
		dm.superModel[0].bc_mask[local_index] = YES;
	}

	//overwrite some of the boundary
	double x_coord, y_coord;
	int id;
	for (int i=0; i<nnodes; i++){
		//mark the boundary only
		x_coord = dm.grid->node[i].x;
		y_coord = dm.grid->node[i].y;
		id = i;
		//id = grid->node[i].id;
		//better, but then permtab always has to be there
		//id = grid->node[grid->permtab[i]].id;
		//printf("x %f, y %f, ID = %d\n",x_coord,y_coord,id);

		if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
			dm.superModel[0].dirichlet_data[id] = 1 + x_coord*x_coord + 2 * y_coord*y_coord;
		}else{
			dm.superModel[0].bc_mask[id]=NO;
		}
	}

	//call fe_newton
	fe_newton(sm); 
	//compare with analytic solution
	//it is a scalar
	double *u_exact;
	u_exact = (double*) tl_alloc(sizeof(double),nnodes);
	compute_exact_solution_poisson(u_exact, nnodes, dm.grid);
//	for(int i=0;i<nnodes;i++){
//		printf("Exact solution[%d] = %f\n",i,u_exact[i]);
//	}
	




	//compute L2 and Linf error
	double l2_err =  l2_error(dm.superModel[0].sol, u_exact, nnodes);
	double linf_err =  linf_error(dm.superModel[0].sol, u_exact, nnodes);

	printf("Final errors: %6.4e , %6.4e\n", l2_err,linf_err);

	//plot grid in h5?
//    strcpy(sm->grid->filename, "residtest");
//    init_hdf5_file(sm->grid);
//    printf("hdf5 initialized\n");
//    sgrid_write_hdf5(sm->grid);
//    printf("hdf5 written\n");
//    sgrid_write_xdmf(sm->grid);
//    printf("xmf written\n");

	//return -1 if failed, 0 if good
	
	if(l2_err < LINEAR_NEWTON_TEST_TOL && linf_err < LINEAR_NEWTON_TEST_TOL){
		err_code=0;
	}
	//printf("Final error code %d\n",err_code);

	//free memory
	u_exact = (double *) tl_free(sizeof(double), nnodes, u_exact);


	smodel_design_free(&dm);


	return err_code;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes analytic solution used to calculate error
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void compute_exact_solution_poisson(double *u_exact, int ndof, SGRID *grid){

	//test case comes from
	//https://jsdokken.com/dolfinx-tutorial/chapter1/fundamentals.html

	//problem is -/\u = f on Omega
	//f = -6
	//u_D = u_exact on dOmega
	//u_exact = 1 + x^2 + 2y^2

	//works for cg only at the moment, would need cell by cell loop for dg
	for(int i =0; i<ndof ; i++){
		u_exact[i] = 1.0 + grid->node[i].x*grid->node[i].x + 2*grid->node[i].y*grid->node[i].y;
	}

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the Newton solver using a nonlinear 
 *   Poisson equation with analytic solution.
 *  equation is -\/.(q(u)\/u) + f = 0
 *                q(u) = 1 + u^2
 * 	              f = 10 + 10x + 20y
 *                Omega = (0,1) x (0,1)
 * 				  u = u_{exact} on \partial \Omega
 *                u_{exact} = 1 + x + 2y 
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int nonlinear_newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax) {
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
	int err_code=-1;
    *(dm.grid) = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );
    int nnodes;
    nnodes = dm.grid->nnodes;
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
    
    //sarray_init_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs);
    sarray_init_value_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs, 1.0);
    dm.superModel[0].LINEAR_PROBLEM = NO;
    dm.superModel[0].max_nonlin_it = 35;
    dm.superModel[0].tol_nonlin = 1e-10;
	dm.superModel[0].inc_nonlin = 1e-8;

 	SMODEL_SUPER *sm;
	sm = &(dm.superModel[0]);

	printf("SETTING UP BCMASK\n");
	double x_coord, y_coord;
	// intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<sm->ndofs; local_index++){
		x_coord = dm.grid->node[local_index].x;
		y_coord = dm.grid->node[local_index].y;
		dm.superModel[0].dirichlet_data[local_index] = 0.0;
		//make initial guess close to true solution
		dm.superModel[0].sol_old[local_index] = 1.0 + x_coord + 2*y_coord + 0.5;
		dm.superModel[0].sol[local_index] = 1.0 + x_coord + 2*y_coord + 0.5;
		dm.superModel[0].lin_sys->dsol[local_index] = 0.0;
		dm.superModel[0].bc_mask[local_index] = YES;
	}

	//overwrite some of the boundary
	
	int id;
	for (int i=0; i<nnodes; i++){
		//mark the boundary only
		x_coord = dm.grid->node[i].x;
		y_coord = dm.grid->node[i].y;
		id = i;
		//id = grid->node[i].id;
		//better, but then permtab always has to be there
		//id = grid->node[grid->permtab[i]].id;
		//printf("x %f, y %f, ID = %d\n",x_coord,y_coord,id);

		if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
			dm.superModel[0].dirichlet_data[id] = 1.0 + x_coord + 2 * y_coord;
		}else{
			dm.superModel[0].bc_mask[id]=NO;
		}
	}
	printf("BCMASK COMPLETE\n");

	//call fe_newton
	fe_newton(sm); 
	//compare with analytic solution
	//it is a scalar
	double *u_exact;
	u_exact = (double*) tl_alloc(sizeof(double),nnodes);
	compute_exact_solution_nonlinear_poisson(u_exact, nnodes, dm.grid);
//	for(int i=0;i<nnodes;i++){
//		printf("Exact solution[%d] = %f\n",i,u_exact[i]);
//	}
	


	//compute L2 and Linf error
	double l2_err =  l2_error(dm.superModel[0].sol, u_exact, nnodes)/nnodes;
	double linf_err =  linf_error(dm.superModel[0].sol, u_exact, nnodes)/nnodes;

	printf("Final errors: %6.4e , %6.4e\n", l2_err,linf_err);


	//return -1 if failed, 0 if good
	
	if(l2_err < NONLINEAR_NEWTON_TEST_TOL && linf_err < NONLINEAR_NEWTON_TEST_TOL){
		err_code=0;
	}
	//printf("Final error code %d\n",err_code);

	//free memory
	u_exact = (double *) tl_free(sizeof(double), nnodes, u_exact);


	smodel_design_free(&dm);


	return err_code;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes analytic solution used to calculate error
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void compute_exact_solution_nonlinear_poisson(double *u_exact, int ndof, SGRID *grid){

	//test case comes from
	//https://jsdokken.com/dolfinx-tutorial/chapter1/fundamentals.html

	//problem is -/\u = f on Omega
	//f = -6
	//u_D = u_exact on dOmega
	//u_exact = 1 + x^2 + 2y^2

	//works for cg only at the moment, would need cell by cell loop for dg
	for(int i =0; i<ndof ; i++){
		u_exact[i] = 1.0 + grid->node[i].x + 2*grid->node[i].y;
	}



}