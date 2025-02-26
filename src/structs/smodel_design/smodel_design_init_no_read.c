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
    int nSuperModels, int nphysics_mats[], char ***elemVarCode, int **coverage_arrays) {

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
        smodel_super_no_read(&(dmod->superModel[imono]), elemVarCode[imono], nphysics_mats[imono], coverage_arrays[imono]);
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
        slin_sys_init_ptrs(&(dmod->lin_sys[i]), &(dmod->superModel[sm_id].my_ndofs),
            &(dmod->superModel[sm_id].ndofs),&(dmod->superModel[sm_id].macro_ndofs), 
            &(dmod->superModel[sm_id].my_ndofs_old), &(dmod->superModel[sm_id].ndofs_old), 
            &(dmod->superModel[sm_id].macro_ndofs_old), 0, dmod->superModel[sm_id].ndofs, 0);

        //NOT IMPLEMENTED, NEED FOR MPI
        //slin_sys_init_ghosts(&(dm->lin_sys[i]), dm->grid, dm->superModel[j].dof_map_local);
        

        //if mono maybe call one routine and if not call another?
        //for now just Mono
        slin_sys_init_sparsity_mono(&(dmod->lin_sys[i]), dmod->superModel[j].elem3d_physics_mat, 
        dmod->superModel[j].elem2d_physics_mat , dmod->superModel[j].elem1d_physics_mat,
        dmod->superModel[j].mat_physics_elem, dmod->grid, dmod->superModel[j].ivars);

//#ifdef _PETSC
//        dm->lin_sys[i].A = PETSC_NULLPTR;
//        dm->lin_sys[i].ksp = PETSC_NULLPTR;
//        dm->lin_sys[i].B = PETSC_NULLPTR;
//        dm->lin_sys[i].X = PETSC_NULLPTR;
//        slin_sys_allocate_petsc_objects(&(dm->lin_sys[i]));
//#endif


    }

//     for(i=0;i<dm->nUnique;i++){
        
//         //define some pointers in each supermodel
//         //assign some hard coded values first, then set the pointers
//         //in general these would require info taken from each superModel

//         //MARK IS SWITCHING FOR DEBUG, PLEASE SWITCH BACK!!!
//         //hard code to only point to first supermodel
//         dm->unique_id[i] = 0; 
//         j= dm->unique_id[i];

//         //allocate the linear systems, need to wrap this into 1 routine
//         // separate into 3 steps for now, how do we want to do ghosts? 
//         //need way to set ghosts before calling sparsity init
//         slin_sys_init_ptrs(&(dm->lin_sys[i]), &(dm->my_ndofs[i]),&(dm->ndofs[i]),&(dm->macro_ndofs[i]),
//         &(dm->my_ndofs_old[i]), &(dm->ndofs_old[i]), &(dm->macro_ndofs_old[i]),
//         0, dm->ndofs[i], 0);
//         //likely will require the grid and maybe fmap? maybe if fmap empty then ghosts is easier
//         slin_sys_init_ghosts(&(dm->lin_sys[i]), dm->grid, dm->superModel[j].dof_map_local);
//         //if mono maybe call one routine and if not call another?
//         //for now just Mono
// //        slin_sys_init_sparsity_mono(&(dm->lin_sys[i]), dm->superModel[j].elem3d_physics_mat_id, 
// //        dm->superModel[j].elem2d_physics_mat_id , dm->superModel[j].elem1d_physics_mat_id ,
// //        dm->superModel[j].elem3d_physics_mat, dm->superModel[j].elem2d_physics_mat,
// //        dm->superModel[j].elem1d_physics_mat ,dm->superModel[j].node_physics_mat_id,
// //        dm->superModel[j].node_physics_mat, dm->grid, dm->superModel[j].dof_map_local);
//         slin_sys_init_sparsity_mono(&(dm->lin_sys[i]), dm->superModel[j].elem3d_physics_mat_id, 
//         dm->superModel[j].elem2d_physics_mat_id , dm->superModel[j].elem1d_physics_mat_id ,
//         dm->superModel[j].elem3d_physics_mat, dm->superModel[j].elem2d_physics_mat,
//         dm->superModel[j].elem1d_physics_mat ,dm->superModel[j].node_physics_mat, 
//         dm->grid, dm->superModel[j].dof_map_local);
// #ifdef _PETSC
//         dm->lin_sys[i].A = PETSC_NULLPTR;
//         dm->lin_sys[i].ksp = PETSC_NULLPTR;
//         dm->lin_sys[i].B = PETSC_NULLPTR;
//         dm->lin_sys[i].X = PETSC_NULLPTR;
//         slin_sys_allocate_petsc_objects(&(dm->lin_sys[i]));
// #endif
//     }

//     //loop through every super model to allocate and assign pointters
//     for(i=0;i<dm->nSuperModels;i++){
//         j = dm->lin_sys_id[i];
//         ndof_temp = dm->ndofs[j];

//         dm->superModel[i].my_ndofs = &(dm->my_ndofs[j]); //pointers to design model, not arrays
//         dm->superModel[i].my_ndofs_old = &(dm->my_ndofs_old[j]);
//         dm->superModel[i].ndofs = &(dm->ndofs[j]);
//         dm->superModel[i].ndofs_old = &(dm->ndofs_old[j]);
//         dm->superModel[i].macro_ndofs = &(dm->macro_ndofs[j]);
//         dm->superModel[i].macro_ndofs_old = &(dm->macro_ndofs_old[j]);
//         dm->superModel[i].bc_mask = (int*) tl_alloc(sizeof(int), ndof_temp);
//         dm->superModel[i].dirichlet_data = (double*) tl_alloc(sizeof(double), ndof_temp);
//         dm->superModel[i].sol = (double*) tl_alloc(sizeof(double), ndof_temp );
//         dm->superModel[i].sol_old = (double*) tl_alloc(sizeof(double), ndof_temp);
//         dm->superModel[i].sol_older = (double*) tl_alloc(sizeof(double), ndof_temp);
//         sarray_init_dbl(dm->superModel[i].sol,ndof_temp);
//         sarray_init_dbl(dm->superModel[i].sol_old,ndof_temp);
//         sarray_init_dbl(dm->superModel[i].sol_older,ndof_temp);

//     }

}
