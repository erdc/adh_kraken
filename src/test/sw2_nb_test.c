/*! \file sw2_nb_test.c This file tests the sw2 engine */
#include "adh.h"
static double NEWTON_TEST_TOL = 2e-3;
static int NEWTON_TEST_NX = 101;//16;
static int NEWTON_TEST_NY = 6;//6;
static double write_testcase_error_nb(SMODEL_SUPER *mod);
static void permute_array(double *arr,int *p, int n);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests a basic sw2 case where we have sloping beach and still
 *  conditions. Mass conservation is tested
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sw2_nb_test(int nt) {
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
	double xmin = 0.0;
	double xmax = 1000.0;
	double ymin = 0.0;
	double ymax = 50.0;
	int npx = NEWTON_TEST_NX;
	int npy = NEWTON_TEST_NY;
	double theta = 45.0;
	double dz = 1.0;
	double a0 = 0.0;
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
    int nnodes;
    nnodes = dm.grid->nnodes;

    //append 1d elements for BC enforcement
    //create a series and string for boundary data without file read
	//first add a string, this will come from smat physics file
	//! Downstream edges
	//EGS 1 2 2 
	//EGS 2 3 2 
	//EGS 3 4 2 
	//EGS 4 5 2 
	//EGS 5 6 2 
	//! Upstream edges
	//EGS 601 602 3 
	//EGS 602 603 3 
	//EGS 603 604 3 
	//EGS 604 605 3 
	//EGS 605 606 3
	dm.grid->nelems1d = 10;
	dm.grid->max_nelems1d = 10;
	dm.grid->macro_nelems1d = dm.grid->nelems1d;
	dm.grid->orig_macro_nelems1d = dm.grid->macro_nelems1d;
	dm.grid->my_nelems1d = dm.grid->nelems1d;
	printf("Attempting to add 1d elements: %d\n",dm.grid->nelems1d);
	dm.grid->elem1d = NULL;
	selem1d_init_alloc_array( &(dm.grid->elem1d), dm.grid->nelems1d);
	printf("Addied 1d elements\n");


	//hard coded, in reality this is given in geo now
	int nodes_segment[2];
	int nnodes_on_elem = 2;
	int n_downstream = 5;
	int n_upstream = 5;
	SVECT nds[nnodes_on_elem]; svect_init_array(nds, nnodes_on_elem);

	//Downstream Edges
	int node_ids[2] = {5,4};
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[0]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 4; node_ids[1] = 3;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[1]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 3; node_ids[1] = 2;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[2]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 2; node_ids[1] = 1;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[3]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 1; node_ids[1] = 0;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[4]), 0, 0, 2, node_ids, 0, nds, 0);

	//Downstream Edges
	node_ids[0] = 600; node_ids[1] = 601;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[5]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 601; node_ids[1] = 602;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[6]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 602; node_ids[1] = 603;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[7]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 603; node_ids[1] = 604;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[8]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 604; node_ids[1] = 605;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = dm.grid->node[node_ids[i]].x;
    	nds[i].y = dm.grid->node[node_ids[i]].y;
    	nds[i].z = dm.grid->node[node_ids[i]].z;
    }
	selem1d_load(&(dm.grid->elem1d[9]), 0, 0, 2, node_ids, 0, nds, 0);


	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Reorder grid to minimize bandwidth
    //++++++++++++++++++++++++++++++++++++++++++++++
    //sgrid_reorder(dm.grid,2);


	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Fill in a simple design model
    //++++++++++++++++++++++++++++++++++++++++++++++
    double tf = 10000.0;
	double dt = tf/ nt;
	double t0 = 0.0;
	
	int nSuperModels = 1;
	int *nphysics_mats;
	nphysics_mats = (int*) tl_alloc(sizeof(int), nSuperModels);
	nphysics_mats[0] = 3;
	//elemVarCode in general is triple pointer, nSuperModels x nphyics_mat[i] x 10
	char ***elemVarCode;
	elemVarCode = (char ***) tl_alloc(sizeof(char**),nSuperModels);
	for (int i = 0; i< nSuperModels; i ++){
		elemVarCode[i] = (char **) tl_alloc(sizeof(char *), nphysics_mats[i]);
		for (int j = 0; j< nphysics_mats[i]; j++){
			elemVarCode[i][j] = (char *) tl_alloc(sizeof(char), 10);
		}
	}
	strcpy(&elemVarCode[0][0][0],"2"); // SW2D
	strcpy(&elemVarCode[0][0][1],"0"); // GW
	strcpy(&elemVarCode[0][0][2],"0"); // Transport

	strcpy(&elemVarCode[0][1][0],"2"); // SW2D
	strcpy(&elemVarCode[0][1][1],"0"); // GW
	strcpy(&elemVarCode[0][1][2],"0"); // Transport

	strcpy(&elemVarCode[0][2][0],"2"); // SW2D
	strcpy(&elemVarCode[0][2][1],"0"); // GW
	strcpy(&elemVarCode[0][2][2],"0"); // Transport


	//mat ids
	int **mat_ids;
	mat_ids = (int **) tl_alloc(sizeof(int*),nSuperModels);
	int nelems = dm.grid->nelems3d + dm.grid->nelems2d + dm.grid->nelems1d;

	for (int i = 0; i < nSuperModels; i++){
		mat_ids[i] = tl_alloc(sizeof(int), nelems);

		//fill in 0 for 2d elem
		//1 for 1d elem
		for (int j = 0; j < dm.grid->nelems2d; j++ ){
			mat_ids[i][j] = 0;
		}
		
		for (int j = dm.grid->nelems2d; j< (dm.grid->nelems2d + n_downstream); j++){
			mat_ids[i][j] = 1;
		}

		for (int j = (dm.grid->nelems2d + n_downstream); j < nelems; j++){
			mat_ids[i][j] = 2;
		}
	}

    smodel_design_init_no_read(&dm, dt, t0, tf, nSuperModels, nphysics_mats, elemVarCode, mat_ids);
    
	
	//use pointer for short hand
	SMODEL_SUPER *sm = &(dm.superModel[0]);
    
	//overwrite inflow
    sm->mat_physics_elem[sm->elem1d_physics_mat[0]].model[0].physics = SW2_BC_FLUX_;
    sm->mat_physics_elem[sm->elem1d_physics_mat[0]].model[0].physics_init = UNSET_INT;

    //overwrite outflow
    sm->mat_physics_elem[sm->elem1d_physics_mat[5]].model[0].physics = SW2_BC_H_;
    sm->mat_physics_elem[sm->elem1d_physics_mat[5]].model[0].physics_init = UNSET_INT;

    //hack together the sw2 structure
    //allocate
    sm->sw = (SSW*) tl_alloc(sizeof(SSW), 1);
    //use an alias
    SSW *sw = dm.superModel[0].sw; //alias for convenience
    //set it up
    printf("allocating ssw\n");
    

    //new way
    //ssw_alloc_init(sw, grid->nnodes, grid->nnodes, grid->nelems2d, 0, 0, 1);
    //Old way
    ssw_alloc_init(sw);
    //set up the dvar map
    printf("allocating sdvar\n");
    sdvar_alloc_init(&(sw->dvar), dm.grid->nnodes, 0, 0, 1, dm.grid->nnodes, dm.grid->nelems2d);
    //hard code set index for wd flag
    sw->WD_FLAG = 0;


    //must also set up the dacont extra terms
    //could be clever hear and pick nnode on each element
    sw->elem_rhs_dacont_extra_terms = (double **) tl_alloc(sizeof(double *), dm.grid->nelems2d);
    for (int i=0; i<dm.grid->nelems2d; i++) {
        sw->elem_rhs_dacont_extra_terms[i] = (double *) tl_alloc(sizeof(double), dm.grid->elem2d[i].nnodes);
        for (int j=0; j<dm.grid->elem2d[i].nnodes; j++) {
        	sw->elem_rhs_dacont_extra_terms[i][j] = 0.;
        }
    }
    //printf("sw vals %d\n",sw->WD_FLAG);
    printf("dvar?? %d\n",sw->dvar.n_dvar_elem_int);
    //initial conditions and things
    // intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index< sm->ndofs; local_index++){
		sm->dirichlet_data[local_index] = 0.0;
		sm->sol_old[local_index] = 0.0;
		sm->sol[local_index] = 0.0;
		sm->lin_sys->dsol[local_index] = 0.0;
		sm->bc_mask[local_index] = NO;
	}

	//overwrite intial condition
	double x_coord, y_coord, z_coord;
	int id;
	for (int i=0; i<nnodes; i++){
		//mark the boundary only
		x_coord = dm.grid->node[i].x;
		y_coord = dm.grid->node[i].y;
		z_coord = dm.grid->node[i].z;
		//id = grid->node[i].id;
		id=i;
		//need to set IC
		sm->sol[id*3] = 1.0 - z_coord;
		sm->sol_old[id*3] = 1.0 - z_coord;
		sm->sol_older[id*3] = 1.0 - z_coord;
		//velocity too
		sm->sol[id*3+1] = 0.070710678;
		sm->sol[id*3+2] = 0.070710678;
		sm->sol_old[id*3+1] = 0.070710678;
		sm->sol_older[id*3+2] = 0.070710678;
		sm->dirichlet_data[id*3] = 1.0 - z_coord;
		sm->dirichlet_data[id*3+1] = 0.070710678;
		sm->dirichlet_data[id*3+2] = 0.070710678;

		//overwrite to 0
//		sm->sol[id*3+1] = 0.0;
//		sm->sol[id*3+2] = 0.0;
//		sm->sol_old[id*3+1] = 0.0;
//		sm->sol_older[id*3+2] = 0.0;
//		sm->dirichlet_data[id*3] = 1.0 - z_coord;
//		sm->dirichlet_data[id*3+1] = 0.0;
//		sm->dirichlet_data[id*3+2] = 0.0;
		//no dirichlet condition
		//if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
		//	continue;
		//}else{
		//	dm.superModel[0].bc_mask[id*3+1]=NO;
		//}
		//printf("Dirichlet data node[%d] = %f\n", i, sm.dirichlet_data[i*3+1]);
	}

	//create a series and string for boundary data without file read
	//first add a string, this will come from smat physics file
	//! Downstream edges
	//EGS 1 2 2 
	//EGS 2 3 2 
	//EGS 3 4 2 
	//EGS 4 5 2 
	//EGS 5 6 2 
	//! Upstream edges
	//EGS 601 602 3 
	//EGS 602 603 3 
	//EGS 603 604 3 
	//EGS 604 605 3 
	//EGS 605 606 3
	//add this to 

	//create a series
	SSERIES *series;
	int nentry = 2; 
	sseries_alloc(&series, nentry, TIME_SERIES, 0);
	series->id = 0;
	series->infact = 1;
	series->outfact = 1;
	series->entry[0].time = 0.0;
  	series->entry[0].time *= series->infact;
  	series->entry[0].value[0] = 0.1;
  	series->entry[1].time = 15000;
    series->entry[1].time *= series->infact;
    series->entry[1].value[0] = 0.1;
    // sanity check the series
  	sseries_check(*series);
  	printf("attempting to add series to linked list\n");
  	sm->series_head = NULL;
  	sm->series_curr = NULL;
  	// add to linked list (which will allocate another, so we can delete this one)
    sseries_add(series, &(sm->series_head), &(sm->series_curr), TRUE); 
	series->id = 1;
	series->infact = 1;
	series->outfact = 1;
	series->entry[0].time = 0.0;
  	series->entry[0].time *= series->infact;
  	series->entry[0].value[0] = 1.0;
  	series->entry[1].time = 15000;
    series->entry[1].time *= series->infact;
    series->entry[1].value[0] = 1.0;
    sseries_add(series, &(sm->series_head), &(sm->series_curr), TRUE); 

    sseries_free(series);

    //create string that points to it!
    //this will be inside elem mat physics now!
    sm->mat_physics_elem[sm->elem1d_physics_mat[0]].bc_ids = (int*) tl_alloc(sizeof(int), 1);
	sm->mat_physics_elem[sm->elem1d_physics_mat[0]].bc_ids[0] = 0;   

	sm->mat_physics_elem[sm->elem1d_physics_mat[5]].bc_ids = (int*) tl_alloc(sizeof(int), 1);
	sm->mat_physics_elem[sm->elem1d_physics_mat[5]].bc_ids[0] = 1;
 

	printf("Initial condition: \n");
	sarray_printScreen_dbl(sm->sol,sm->ndofs, "sol");


	printf("Calling time loop\n");
	//set forward_step and call timeloop
	time_loop(&dm); 

	//get initial mass
	double *init_head;
    init_head = (double *) tl_alloc(sizeof(double), nnodes);
    SVECT2D *init_vel;
    init_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    //fill in temp arrays
    for (int inode=0; inode<nnodes; inode++){
    	init_head[inode] = dm.superModel[0].dirichlet_data[inode*3];
    	init_vel[inode].x = dm.superModel[0].dirichlet_data[inode*3+1];
    	init_vel[inode].y = dm.superModel[0].dirichlet_data[inode*3+2];
    }

	//printf("Initial grid mass = %f\n",initial_grid_mass);
	//extract second variable here
	double total_error = write_testcase_error_nb(&(dm.superModel[0])); 
	printf("Final error: %6.4e\n", total_error);
	//plot grid in h5?
//    strcpy(sm.grid->filename, "residtest");
//    init_hdf5_file(sm.grid);
//    printf("hdf5 initialized\n");
//    sgrid_write_hdf5(sm.grid);
//    printf("hdf5 written\n");
//    sgrid_write_xdmf(sm.grid);
//    printf("xmf written\n");

	//return -1 if failed, 0 if good
	
	if(total_error < NEWTON_TEST_TOL){
		ierr=0;
	}
	//printf("Final error code %d\n",err_code);




	smodel_design_free(&dm);


	return ierr;
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
// note :: may get natural errors from water head not remaining flat at in flow edge
double write_testcase_error_nb(SMODEL_SUPER *mod) {
    
    int inode;
    double convert_to_rads = 3.141592653589793 / 180.;
    
    // the inflow is 0.1 m/s at a 45 degree angle
    double analytic_vx = 0.1/sqrt(2.);
    double analytic_vy = analytic_vx;
    
    double error_vx = 0., error_vx_max;
    double error_vy = 0., error_vy_max;
    double error_h = 0., error_h_max;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sol[inode*3+1] - analytic_vx);
        error_vy += fabs(mod->sol[inode*3+2] - analytic_vy);
        error_h  += fabs(mod->sol[inode*3] - 1.);  // this is just to see how far from rigid lid
    }

    return error_vx/mod->grid->macro_nnodes/analytic_vx + error_vy/mod->grid->macro_nnodes/analytic_vy + error_h/mod->grid->macro_nnodes/1.;
}