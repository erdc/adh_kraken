/*! \file residual_test.c This file tests the residual assembly */
#include "adh.h"
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
int residual_test(void) {
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
	double xmin = 0.0;
	double xmax = 2.0;
	double ymin = 0.0;
	double ymax = 2.0;
	int npx = 13;
	int npy = 5;
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
	strcpy(&elemVarCode[0][0][0],"2");//SW2D
	strcpy(&elemVarCode[0][0][1],"0"); //GW
	strcpy(&elemVarCode[0][0][2],"0"); //Transport

	printf("elem var code %s\n",elemVarCode[0][0]);
	//mat ids
	int **mat_ids;
	mat_ids = (int **) tl_alloc(sizeof(int*),nSuperModels);
	int nelems = dm.grid->nelems2d + dm.grid->nelems1d + dm.grid->nelems1d;
	for (int i = 0; i < nSuperModels; i++){
		mat_ids[i] = tl_alloc(sizeof(int), nelems);
		sarray_init_int(mat_ids[i],nelems);
	}

    smodel_design_init_no_read(&dm, dt, t0, tf, nSuperModels, nphysics_mats, elemVarCode, mat_ids);
    //smodel_design_no_read_simple(&dm, dt, t0, tf, 1, elemVarCode, grid);
	////specify elemental physics and other properties in super model


//	char elemVarCode[4]; 
//	strcpy(&elemVarCode[0],"2");//SW2D
//	strcpy(&elemVarCode[1],"0"); //GW
//	strcpy(&elemVarCode[2],"0"); //Transport
//	//printf("GRID NELEMS2D = %d\n",grid.nelems2d);
//	//smodel_super_no_read_simple(&sm, dt, t0, tf, 0 , 1, 0, elemVarCode);
//	smodel_design_no_read_simple(&dm, dt, t0, tf, 1, elemVarCode, grid);
//	//printf("NDOFS %d\n",dm->ndofs[0]);//

////	//assemble a residual and check correctness
//	assemble_residual(&(dm.superModel[0]), dm.grid);

	//print final residual
//	for(int local_index=0;local_index<grid.nnodes;local_index++){
//		printf("Node %d: (x,y) = {%f,%f}, Residual = {%f,%f,%f}\n",grid.node[local_index].gid,grid.node[local_index].x,grid.node[local_index].y,sm.residual[local_index*3],sm.residual[local_index*3+1],sm.residual[local_index*3+2]);
//	}




    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Free Memory
    //++++++++++++++++++++++++++++++++++++++++++++++
    smodel_design_free(&dm);
    
	

	return 0;
}
