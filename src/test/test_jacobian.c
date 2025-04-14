/*! \file test_jacobian.c This file tests the assembly of a Jacobian matrix */
#include "adh.h"
static int check_jacobian(SMODEL_SUPER *sm);
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
int test_jacobian(int npx, int npy, double xmin, double xmax, double ymin, double ymax) {
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
    sgrid_create_node_to_node_graph(dm.grid);
    sgrid_reorder(dm.grid,2);
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Fill in a simple design model
    //++++++++++++++++++++++++++++++++++++++++++++++
    double dt = 1.0;
	double t0 = 0.0;
	double tf = 1.0;
	int nSuperModels = 1;
	int *nphysics_mats = (int*) tl_alloc(sizeof(int), nSuperModels);
	nphysics_mats[0] = 1;

	//need to allocate the npdes argument, this will be a double pointer
	// which will be nSuperModels x nphysics_mat[i]
	// each entry is number of pdes per physics material
	int **npdes = (int**) tl_alloc(sizeof(int*), nSuperModels);
	int *npdes_total = (int *) tl_alloc(sizeof(int), nSuperModels);
	sarray_init_int(npdes_total, nSuperModels);
	for(int i = 0 ; i < nSuperModels ; i++){
		npdes[i] = (int*) tl_alloc(sizeof(int), nphysics_mats[i]);
		for(int j = 0; j < nphysics_mats[i]; j++){
			//this case just one pde on each material within each super model
			npdes[i][j] = 1;
			npdes_total[i] += npdes[i][j];
		}
	}

	//use this to create modelvsbc double pointer
	//this is an array that is nSuperModels x npdes_total[i]
	//where npdes_total[i] = sum_j(npdes[i][j])
	int **modelvsbc = (int**) tl_alloc(sizeof(int*), nSuperModels);
	int *nmodel = (int*) tl_alloc(sizeof(int), nSuperModels);
	sarray_init_int(nmodel, nSuperModels);
	for (int i = 0; i < nSuperModels; i++){
		modelvsbc[i] = tl_alloc(sizeof(int), npdes_total[i]);
		for (int j = 0 ; j < npdes_total[i]; j++){
			//1 for model, 0 for bc
			//this test no bc integrals, only models
			modelvsbc[i][j] = 1;
			nmodel[i]+=modelvsbc[i][j];
		}
	}

	//Now save strings to be read, this will determine pdes to be active
	//first do model strings
	//model_strings in general is triple pointer, nSuperModels x nmodel[i] x max_char
	// where nmodel[i] = sum_j(modelvsbc[i][j])
	// the models are assumed to be in order that materials are ordered
	// so if mat 1 has 2 pdes and mat 2 has 1 pde
	// model_strings[0] and model_strings[1] will be assigned to mat 1
	// and model_strings[2] will be assigned to mat 2

	int max_char = 30;
	char ***model_strings = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	char ***bc_phystype = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	char ***bc_type = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	char ***bc_vartype = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	int **bc_iseries = (int **) tl_alloc(sizeof(int*), nSuperModels);

	for (int i = 0; i< nSuperModels; i ++){
		model_strings[i] = (char **) tl_alloc(sizeof(char *), nmodel[i]);
		bc_phystype[i] = NULL;
		bc_type[i] = NULL;
		bc_vartype[i] = NULL;
		bc_iseries[i] = NULL;
		for (int j = 0; j< nmodel[j]; j++){
			model_strings[i][j] = (char *) tl_alloc(sizeof(char), max_char);
		}
	}
	//this particular test case is only 1 super model and 1 material that is POISSON2D
	strcpy(model_strings[0][0],"POISSON2D"); // POISSON


	//mat ids
	int **mat_ids;
	mat_ids = (int **) tl_alloc(sizeof(int*),nSuperModels);
	int nelems = dm.grid->nelems3d + dm.grid->nelems2d + dm.grid->nelems1d;
	for (int i = 0; i < nSuperModels; i++){
		mat_ids[i] = tl_alloc(sizeof(int), nelems);
		sarray_init_int(mat_ids[i],nelems);
	}

    
    smodel_design_init_no_read(&dm, dt, t0, tf, nSuperModels, nphysics_mats, npdes, modelvsbc,
    	model_strings, bc_phystype, bc_type, bc_vartype, bc_iseries, mat_ids);
    
    //sarray_init_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs);
    sarray_init_value_dbl(dm.superModel[0].sol, dm.superModel[0].ndofs, 1.0);
    dm.superModel[0].LINEAR_PROBLEM = YES;


	assemble_jacobian(&(dm.superModel[0]));


	ierr = check_jacobian(&(dm.superModel[0]));

	//free stuff
    smodel_design_free(&dm);
	return ierr;
}


int check_jacobian(SMODEL_SUPER *sm){
	int *indptr = sm->lin_sys->indptr_diag;
	int *cols = sm->lin_sys->cols_diag;
	double *vals = sm->lin_sys->vals_diag;
	int nrow = *(sm->lin_sys->local_size);
	int nnz = sm->lin_sys->nnz_diag;
	int n_connections;
	int NNZ = 0;
	int err = 0;
	//check if number of columns on each row makes sense
	int nnodes = sm->grid->nnodes;
	int nd1, nd2, nd3;
	int *nodes;
	int *n_con;
	int nnodes_on_elem;
	int current_node;
	int other_node;
	n_con = (int *) tl_alloc(sizeof(int), nnodes);
	sarray_init_int(n_con, nnodes);

	int i,j,k;

	//we know max nodes this can be connected to will be 9
	int max_con = 9;
	int **temp_edgetab;
	temp_edgetab = (int**) tl_alloc(sizeof(int*), nnodes);
    for(int j=0;j<nnodes;j++){
        temp_edgetab[j] = (int*) tl_alloc(sizeof(int), max_con);
        for(int k=0;k<max_con;k++){
            temp_edgetab[j][k]=INT_MAX;
        }
    }


    //First set of loops is solely to establish how many nodes are connected to each node
    for (i=0;i<sm->grid->nelems2d;i++){
    	//nnodes on the element
        nnodes_on_elem = sm->grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = sm->grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node and add
            for (k=0;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
      
        }
    }


	for (i=0;i<nnodes;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_edgetab[i], n_con[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        n_connections = sarray_unique_int(temp_edgetab[i], n_con[i]);
        //overwrite nnz row with sarray_unique_int?
        n_con[i] = n_connections;
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        NNZ+=n_connections;
    }



	for (int i = 0;i<nrow;i++){
		if ( (indptr[i+1] - indptr[i]) != n_con[i]){
			err+=1;
		}

	}

	if (err!=0){err=1;} 

	return err;

}