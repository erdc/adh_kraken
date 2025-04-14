/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_init_no_read.c This file collects methods of the SMODEL_DESIGN structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initialize design model without reading a file for testing
 *             Assumes the design grid is already set up
 *             Also assumes we are only dealing with one supermodel (for now)
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod                (SMODEL_SUPER *)  an AdH superModel
 * @param[in]  FILE                    (FILE *) the SuperModel input file
 * \note This supermodel is already assumed to have a grid pointer within it that is populated
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_init_no_read(SMODEL_DESIGN *dmod, double dt_in, double t_init, double t_final,
    int nSuperModels, int nphysics_mats[], int **npdes, int **modelvsbc, char ***model_strings, 
    char ***bc_phystype, char ***bc_type, char ***bc_vartype, int **bc_iseries, int **coverage_arrays) {

    int i,j;
    int isSimple=0;
    int ndof_temp;
    
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Designer Model No Read\n");
        printf("------------------------------------------------------\n");
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Set Time Step Information
    // Simplified set up assuming SI units
    // and fixed time step
    //++++++++++++++++++++++++++++++++++++++++++++++
    dmod->ts.t_init = t_init;
    dmod->ts.t_prev = t_init;
    dmod->ts.t_final = t_final;
    //instead of
    //dm->series_dt = sseries_read_allocate(NULL,NULL,DT_SERIES,&token,UNSET_INT); 
    //do manually
    int nentries = 2;
    double time[2] = {t_init,t_final};
    double values[2] = {dt_in, dt_in};
    sseries_alloc(&(dmod->series_dt), nentries, DT_SERIES, UNSET_INT);
    dmod->series_dt->id = 0;
    dmod->series_dt->nnodes = UNSET_INT;
    dmod->series_dt->size = nentries;
    for (int ientry=0; ientry < dmod->series_dt->size; ientry++) {
        dmod->series_dt->entry[ientry].time  = time[ientry];
        for (int ivalue=0; ivalue < dmod->series_dt->nvalues; ivalue++) {
            dmod->series_dt->entry[ientry].value[ivalue] = values[ientry];
        } 
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Now count SuperModels
    // Provided as input
    //++++++++++++++++++++++++++++++++++++++++++++++
    dmod->nSuperModels = nSuperModels;


    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate and Initialize SuperModels
    //++++++++++++++++++++++++++++++++++++++++++++++
    smodel_super_alloc_init(&dmod->superModel,dmod->nSuperModels);

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Assign designer grid to each superModel
    //++++++++++++++++++++++++++++++++++++++++++++++
    for (int imono=0; imono<dmod->nSuperModels; imono++) {
        dmod->superModel[imono].grid = dmod->grid;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // set superModel parameters
    //++++++++++++++++++++++++++++++++++++++++++++++
    for (int imono=0; imono<dmod->nSuperModels; imono++) {
        if (DEBUG) {
            printf("------------------------------------------------------\n");
            printf("------------------------------------------------------\n");
            printf("Initiating Super Model Without File Read\n");
            printf("------------------------------------------------------\n");
        }
        //with read
        //smodel_super_read(&(dmod->superModel[imono]));
        //without read replace with
        printf("attempting to call super model no read once\n");
        smodel_super_no_read(&(dmod->superModel[imono]), nphysics_mats[imono], npdes[imono], modelvsbc[imono],
        model_strings[imono], bc_phystype[imono], bc_type[imono], bc_vartype[imono], bc_iseries[imono],
        coverage_arrays[imono]);
        //how to do this properly?
        dmod->superModel[imono].ts = &(dmod->ts);
        //exit(-1);

//        if (DEBUG) {
//            printf("------------------------------------------------------\n");
//            printf("Reading SuperModel '%s' Initialization File\n",dmod->superModel[imono].filebase);
//            printf("------------------------------------------------------\n");
//        }
        //for now initial conditions must be set outside of this
        //smodel_super_read_init(&(dmod->superModel[imono]),dmod->superModel[imono].filebase);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Open Parameter Coverage File
    // One coverage can be used for multiple superModels
    // For now left unimplemented
    //++++++++++++++++++++++++++++++++++++++++++++++filebase,NULL,NULL,suffix,mode,TRUE);
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Designer Model Parameter File Read\n");
        printf("------------------------------------------------------\n");
    }
    //scoverage_read(&dmod->params,filebase);
    //scoverage_no_read(...,...);


    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate linear systems
    // for now, each superModel has its own
    // linear system, need to change in future
    // for models to share
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Lin Sys\n");
        printf("------------------------------------------------------\n");
    }
    dmod->nlin_sys = dmod->nSuperModels;
    slin_sys_alloc_array(&(dmod->lin_sys), dmod->nlin_sys);
    //in general this will be an nlin_sys array pointing
    //to first supermodel with that sparsity pattern
    dmod->unique_id = (int *) tl_alloc(sizeof(int), dmod->nlin_sys);
    //need routine to set this generally but for now
    sarray_init_range_int(dmod->unique_id, dmod->nlin_sys);
    int sm_id; 
    for (int imono=0; imono<dmod->nlin_sys; imono++){
        sm_id = dmod->unique_id[imono];
        //allocate the linear systems, need to wrap this into 1 routine
        // separate into 3 steps for now, how do we want to do ghosts? 
        //need way to set ghosts BEFORE calling sparsity init

        //WARNING: ONLY SET FOR SERIAL (SEE LAST THREE ARGUMENTS)
        slin_sys_init_ptrs(&(dmod->lin_sys[imono]), &(dmod->superModel[sm_id].my_ndofs),
            &(dmod->superModel[sm_id].ndofs),&(dmod->superModel[sm_id].macro_ndofs), 
            &(dmod->superModel[sm_id].my_ndofs_old), &(dmod->superModel[sm_id].ndofs_old), 
            &(dmod->superModel[sm_id].macro_ndofs_old), 0, dmod->superModel[sm_id].ndofs, 0);

        //NOT IMPLEMENTED, NEED FOR MPI
        //slin_sys_init_ghosts(&(dm->lin_sys[i]), dm->grid, dm->superModel[j].dof_map_local);
        

        //if mono maybe call one routine and if not call another?
        //for now just Mono
        slin_sys_init_sparsity_mono(&(dmod->lin_sys[imono]), dmod->superModel[sm_id].elem3d_physics_mat, 
        dmod->superModel[sm_id].elem2d_physics_mat , dmod->superModel[sm_id].elem1d_physics_mat,
        dmod->superModel[sm_id].mat_physics_elem, dmod->grid, dmod->superModel[sm_id].ivars);

#ifdef _PETSC
        dmod->lin_sys[imono].A = PETSC_NULLPTR;
        dmod->lin_sys[imono].ksp = PETSC_NULLPTR;
        dmod->lin_sys[imono].B = PETSC_NULLPTR;
        dmod->lin_sys[imono].X = PETSC_NULLPTR;
        slin_sys_allocate_petsc_objects(&(dmod->lin_sys[imono]));
#endif


    }
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Associate the correct supermodel to correct
    // linear system
    //++++++++++++++++++++++++++++++++++++++++++++++
    //Ask corey
    dmod->superModel[0].lin_sys = &(dmod->lin_sys[0]);
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Design model init completed\n");
        printf("------------------------------------------------------\n");
    }

}
