/*! \file sw2_wd_test.c This file tests the sw2 engine */
#include "adh.h"
static double NEWTON_TEST_TOL = 1e-7;
static int NEWTON_TEST_NX = 16;
static int NEWTON_TEST_NY = 6;
static double write_testcase_error_wet_dry(SMODEL_SUPER *mod, double initial_grid_mass);
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
int sw2_wd_test(int npx, int npy, int nt) {

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
	double ymax = 500.0;
	double theta = 0.0;
	double dz = 1.0;
	double a0 = 0.0;
	double ax = 0.0015;
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

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Reorder grid to minimize bandwidth
    //++++++++++++++++++++++++++++++++++++++++++++++
    sgrid_reorder(dm.grid,2);

	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Fill in a simple design model
    //++++++++++++++++++++++++++++++++++++++++++++++
    double tf = 864000.0;
	double dt = tf/ nt;
	double t0 = 0.0;
	
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
	strcpy(&elemVarCode[0][0][0],"2"); // SW2D
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
    


    //hack together the sw2 structure
    //allocate
    dm.superModel[0].sw = (SSW*) tl_alloc(sizeof(SSW), 1);
    //use an alias
    SSW *sw = dm.superModel[0].sw; //alias for convenience
    //set it up
    printf("allocating ssw\n");
    ssw_alloc_init(sw);

    //set up the dvar map
    printf("allocating sdvar\n");
    sdvar_alloc_init(&(sw->dvar), nnodes, 0, 0, 1, nnodes, dm.grid->nelems2d);
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
    printf("sw vals %d\n",sw->WD_FLAG);
    printf("dvar?? %d\n",sw->dvar.n_dvar_elem_int);
    //initial conditions and things
    // intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<dm.superModel[0].ndofs; local_index++){
		dm.superModel[0].dirichlet_data[local_index] = 0.0;
		dm.superModel[0].sol_old[local_index] = 0.0;
		dm.superModel[0].sol[local_index] = 0.0;
		dm.superModel[0].lin_sys->dsol[local_index] = 0.0;
		dm.superModel[0].bc_mask[local_index] = NO;
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
		dm.superModel[0].sol[id*3] = 1.0 - z_coord;
		dm.superModel[0].sol_old[id*3] = 1.0 - z_coord;
		dm.superModel[0].sol_older[id*3] = 1.0 - z_coord;
		dm.superModel[0].dirichlet_data[id*3] = 1.0 - z_coord;
		//no dirichlet condition
		//if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
		//	continue;
		//}else{
		//	dm.superModel[0].bc_mask[id*3+1]=NO;
		//}
		//printf("Dirichlet data node[%d] = %f\n", i, sm.dirichlet_data[i*3+1]);
	}

	printf("Calling time loop\n");

	int matid = dm.superModel[0].elem2d_physics_mat[0];
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
    SFLAGS dummy;
	double initial_grid_mass = tl_find_grid_mass_elem2d(dm.superModel[0].density, NULL, NULL, init_head,dm.superModel[0].grid, dummy);
	printf("Initial grid mass = %f\n",initial_grid_mass);

	//extract second variable here
	double total_error = write_testcase_error_wet_dry(&(dm.superModel[0]),initial_grid_mass); 



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
double write_testcase_error_wet_dry(SMODEL_SUPER *mod, double initial_grid_mass) {
    int inode;
    double error_vx = 0., error_vy = 0., error_h = 0., max_head = 0., max_u = 0., max_v = 0., total_time_mass_flux= 0.0;
    double error_vx_max, error_vy_max, error_h_max, max_head_grid, max_u_grid, max_v_grid;

    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sol[inode*3+1]); // water initially at rest
        error_vy += fabs(mod->sol[inode*3+2]); // water initially at rest
        error_h  += fabs(mod->sol[inode*3] - mod->dirichlet_data[inode*3]);
        if (fabs(mod->dirichlet_data[inode*3]) > max_head) max_head = fabs(mod->dirichlet_data[inode*3]);
    }
    //NEED TO DO
    double *head;
    head = (double *) tl_alloc(sizeof(double), mod->grid->nnodes);
    SVECT2D *vel;
    vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), mod->grid->nnodes);
    //fill in temp arrays
    for (inode=0; inode<mod->grid->my_nnodes; inode++){
    	head[inode] = mod->sol[inode*3];
    	vel[inode].x = mod->sol[inode*3+1];
    	vel[inode].y = mod->sol[inode*3+2];
    }
    SFLAGS dummy;
    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, head, vel, mod->grid, dummy, initial_grid_mass, NULL, NULL, mod->ts->dt, &total_time_mass_flux);
    
    printf("Grid mass error %6.4e\n", grid_mass_error);


    head = tl_free(sizeof(double),mod->grid->nnodes,head);
    vel = tl_free(sizeof(SVECT2D),mod->grid->nnodes,vel);

    return error_vx + error_vy + error_h;

//    if (max_head < 1e-6) max_head = 1.;
//    if (max_u < 1e-6) max_u = 1.;
//    if (max_v < 1e-6) max_v = 1.;
//#ifdef _MESSG
//    error_vx_max=messg_dsum(error_vx,mod->grid->smpi->ADH_COMM);
//    error_vy_max=messg_dsum(error_vy,mod->grid->smpi->ADH_COMM);
//    error_h_max=messg_dsum(error_h,mod->grid->smpi->ADH_COMM);
//    max_head_grid=messg_dmax(max_head,mod->grid->smpi->ADH_COMM);
//    max_u_grid=messg_dmax(max_u,mod->grid->smpi->ADH_COMM);
//    max_v_grid=messg_dmax(max_v,mod->grid->smpi->ADH_COMM);
//    error_vx=error_vx_max;
//    error_vy=error_vy_max;
//    error_h=error_h_max;
//    max_head=max_head_grid;
//    max_u=max_u_grid;
//    max_v=max_v_grid;
//#endif
//    if(mod->grid->smpi->myid==0){
//        FILE *fp;
//        fp = fopen("error.out", "w");
//        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/mod->grid->macro_nnodes,error_vx/mod->grid->macro_nnodes/max_u);
//        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->macro_nnodes,error_vy/mod->grid->macro_nnodes/max_v);
//        fprintf(fp,"head abs error: %30.20e :: relative error: %30.20e \n",error_h/mod->grid->macro_nnodes,error_h/mod->grid->macro_nnodes/max_head);
//        fprintf(fp,"grid_mass_error: %30.20e :: relative error: %30.20e \n",grid_mass_error, grid_mass_error/mod->initial_grid_mass);
//        fclose(fp);
//    }
}
void permute_array(double *arr,int *p, int n){
	double *temp;
	temp = (double *) tl_alloc(sizeof(double),n);
	for(int i =0;i<n;i++){
		temp[p[i]] = arr[i];
	}
	// Copy permuted elements back to the original array
    for (int i = 0; i < n; i++) {
        arr[i] = temp[i];
    }
	temp = (double *) tl_free(sizeof(double),n,temp);
}
