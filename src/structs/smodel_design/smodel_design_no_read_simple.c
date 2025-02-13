/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_design_no_read_simple.c This file collects methods of the SMODEL_DESIGN structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initialize super model without reading a file for testing
 *             Only designed to take one material for an entire mesh
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
void smodel_design_no_read_simple(SMODEL_DESIGN *dm, double dt_in, double t_init, double t_final,
    int nphysics_mat, char elemVarCode[4] ,
    SGRID *grid) {
    
    int i,j;
    int isSimple=0;
    int ndof_temp;
    
    printf("Initializing design model without file read\n");
    int nSuperModels = 1;
    smodel_design_alloc_init(&dm,nSuperModels);
    dm->nMono = 1;
    dm->nUnique = 1;
    dm->grid = grid;
    
    dm->lin_sys = (SMODEL_SUPER*) tl_alloc(sizeof(SMODEL_SUPER), dm->nSuperModels);

    dm->ts.dt = dt_in;
    dm->ts.dt_old = dt_in;
    dm->ts.dt_err = dt_in;
    dm->ts.dt_prev = dt_in;
    dm->ts.t_init = t_init;
    dm->ts.t_prev = t_init;
    dm->ts.t_final = t_final;
    dm->ts.t_adpt_flag = 0;

    for(i=0;i<dm->nSuperModels;i++){
        
        
        
        
        
        smodel_super_no_read_simple(&(dm->superModel[i]), &(dm->ts.dt), &(dm->ts.t_init), &(dm->ts.t_prev),
        &(dm->ts.t_final),nphysics_mat, elemVarCode, isSimple,dm->grid, &(dm->lin_sys[i]));
        
        dm->superModel[i].tol_nonlin = 1e-5;
        dm->superModel[i].inc_nonlin = 1e-3;
        dm->superModel[i].max_nonlin_linesearch_cuts = 5;
        dm->superModel[i].it_count_nonlin_failed = 0;
        dm->superModel[i].max_nonlin_it = 25;
        dm->superModel[i].LINEAR_PROBLEM = YES;
        dm->superModel[i].force_nonlin_it = NO;
        dm->superModel[i].force_nonlin_it = NO;
        dm->superModel[i].nonlinear_it_total = 0;
        //hard code to only point to first lin sys
        dm->lin_sys_id[i] = 0; 


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
