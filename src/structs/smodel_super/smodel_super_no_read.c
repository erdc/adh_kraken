/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super_no_read.c This file sets data for a SUPER_MODEL structure  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     sets SuperModel data without file read
 *             used for unit testing
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] sm  (SUPER_MODEL *)  pointer to an AdH superModel
 * @param[in] codes  (char [][10])  array of codes, should have 1 for each mat physics
 * @param[in] nmat_physics (int) the number of physics materials
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_no_read(SMODEL_SUPER *sm, int nmat_physics, int *npdes, int *modelvsbc,
    char **model_strings, char **bc_phystype, char **bc_type, char **bc_vartype, int *bc_iseries,
    int *mat_ids) {
    int i,j;
    SGRID *grid = sm->grid; // alias
    SMAT_PHYSICS *mat = NULL;
    SPDE *pde = NULL;
    int npdes_total = 0;
    char *token;
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Set Element Physics Materials without BC input
    // This section replaces what happens in read_bc_MODEL(sm,fp);
    //++++++++++++++++++++++++++++++++++++++++++++++
    sm->nmat_physics = nmat_physics;
    smat_physics_alloc_init_array(&(sm->mat_physics_elem),sm->nmat_physics);
    for (i = 0; i < sm->nmat_physics; i++ ){
        mat = &(sm->mat_physics_elem[i]);
        mat->id = i;
        mat->npdes = npdes[i];
        mat->pde = (SPDE *) tl_alloc(sizeof(SPDE), mat->npdes);
        for (j=0; j<mat->npdes; j++) {
            mat->pde[j].imod = UNSET_INT;
            mat->pde[j].iseries = UNSET_INT;
            mat->pde[j].resid = NULL;
            mat->pde[j].init  = NULL;
            npdes_total++;
        }
    }

    //fill in SPDE stuff
    int ntrns_mod = 0, ntrns_bc = 0, series_id = UNSET_INT;
    int mat_id = 0;
    int pde_ctr = 0;
    int model_ctr = 0;
    int bc_ctr = 0;
    int itmp;
    mat = &(sm->mat_physics_elem[0]);
    for (i = 0; i <npdes_total; i++){
        //read modelvsbc to know which type of string to get
        //model is 1, bc is 0
        if (modelvsbc[i] == 1){
            model_flag_switch(model_strings[model_ctr],sm->flags.model);
            pde = &mat->pde[pde_ctr];
            pde->resid = select_resid_func(model_strings[model_ctr],&(pde->imod),ntrns_mod);
            if (strcmp(model_strings[model_ctr], "TRNS") == 0) {ntrns_mod++;}
            if (DEBUG) {printf("MODEL found: %s on material: %d || pde->imod: %d\n",model_strings[model_ctr],mat_id,pde->imod);}
            pde_ctr++;
            model_ctr++;
        }else if(modelvsbc[i] == 0){
            pde = &mat->pde[bc_ctr];
            if (strcmp(bc_phystype, "TRNS") == 0) {ntrns_bc++;}
            pde->resid = select_bresid_func(bc_phystype[bc_ctr],bc_type[bc_ctr],
                bc_vartype[bc_ctr],&(pde->imod),ntrns_bc);
            pde->iseries = bc_iseries[bc_ctr];
            if (DEBUG) {printf("BC found: %s on material: %d || pde->imod: %d\n",bc_phystype[bc_ctr],mat_id,pde->imod);}
            pde_ctr++;
            bc_ctr++;
        }else{
            tl_error("modelvsbc must be 1 or 0 (1 for model 0 for bc)");
        }
        //update indices if we have filled out all pdes on this mat
        if (pde_ctr == mat->npdes){
            mat_id++;
            mat = &(sm->mat_physics_elem[mat_id]);
            pde_ctr=0;
        }
    }
    assert(mat_id == sm->nmat_physics);
    smat_physics_set_ivar_pos(sm->mat_physics_elem, sm->nmat_physics);
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Elemental Physics Coverage File
    // Set with array
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (grid->nelems3d > 0) {sm->elem3d_physics_mat = (int *) tl_alloc(sizeof(int),grid->nelems3d);}
    if (grid->nelems2d > 0) {sm->elem2d_physics_mat = (int *) tl_alloc(sizeof(int),grid->nelems2d);}
    if (grid->nelems1d > 0) {sm->elem1d_physics_mat = (int *) tl_alloc(sizeof(int),grid->nelems1d);}
    int n1d=0,n2d=0,n3d=0,nmat1d[30],nmat2d[30],nmat3d[30];
    sarray_init_int(nmat1d,30); 
    sarray_init_int(nmat2d,30); 
    sarray_init_int(nmat3d,30);
    for (i=0;i<grid->nelems3d;i++){
        sm->elem3d_physics_mat[n3d] = mat_ids[i];
        nmat3d[sm->elem3d_physics_mat[n3d]]++;
        n3d++;
    }
    for (i=grid->nelems3d; i<(grid->nelems3d + grid->nelems2d);i++){
        sm->elem2d_physics_mat[n2d] = mat_ids[i];
        nmat2d[sm->elem2d_physics_mat[n2d]]++;
        n2d++;
    }
    for (i = grid->nelems3d + grid->nelems2d; i<(grid->nelems3d + grid->nelems2d + grid->nelems1d);i++){
        sm->elem1d_physics_mat[n1d] = mat_ids[i];
        nmat1d[sm->elem1d_physics_mat[n1d]]++;
        n1d++;
    }
    assert(n1d == sm->grid->nelems1d); 
    assert(n2d == sm->grid->nelems2d);
    assert(n3d == sm->grid->nelems3d);
    if (DEBUG) {
        printf("Physics Coverage File Info:\n");
        printf("-- n1d: %d || nelems1d: %d\n",n1d,sm->grid->nelems1d);
        printf("-- n2d: %d || nelems2d: %d\n",n2d,sm->grid->nelems2d);
        printf("-- n3d: %d || nelems3d: %d\n",n3d,sm->grid->nelems3d);
        for (i=0; i<30; i++) {if (nmat3d[i] > 0) printf("-- 3D physics material[%d] used %d times\n",i,nmat3d[i]);}
        for (i=0; i<30; i++) {if (nmat2d[i] > 0) printf("-- 2D physics material[%d] used %d times\n",i,nmat2d[i]);}
        for (i=0; i<30; i++) {if (nmat1d[i] > 0) printf("-- 1D physics material[%d] used %d times\n",i,nmat1d[i]);}
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Create Node Materials Based by Pointing to Highest DOF SMAT_PHYSICS
    //++++++++++++++++++++++++++++++++++++++++++++++

    if (DEBUG) {
        printf(">creating mat_physics_node\n");
    }
    sm->mat_physics_node = (SMAT_PHYSICS **) tl_alloc(sizeof(SMAT_PHYSICS *),grid->nnodes);
    smat_physics_set_nodal_pointers(sm->mat_physics_node, grid, sm->elem1d_physics_mat,
    sm->elem2d_physics_mat, sm->elem3d_physics_mat, sm->mat_physics_elem);

    
     //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Find all independent variables used in the SuperModel
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">creating ivar_pos\n");
    }
    int FLAGS[adh_def.n_ivars];
    sivar_position_init(&(sm->ivar_pos));
    smat_physics_position_flag(sm->mat_physics_node,grid->nnodes,FLAGS); 
    sivar_position_map(&(sm->ivar_pos),FLAGS);
    //printf("sm->ivar_pos.n: %d\n",sm->ivar_pos.n);
    //printf("sm->ivar_pos.ntrns: %d\n",sm->ivar_pos.ntrns);
    if (DEBUG) {
        printf("Variable Position Info:\n");
        sivar_position_printScreen(&sm->ivar_pos);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate the map from the independent variables 
    // to their position in the solution vector
    // using the sivar position that has been filled out
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">allocating solution variables\n");
    }
    SIVAR_POSITION *ip;
    ip = (SIVAR_POSITION *) tl_alloc(sizeof(SIVAR_POSITION), grid->nnodes);
    for (i=0; i<grid->nnodes; i++) {ip[i] = sm->mat_physics_node[i]->ivar_pos;}
    sm->ndofs = sivar_position_build_dof_map(&sm->ivar_pos,grid->nnodes,ip,&sm->ivars);
    if (DEBUG) {
        printf(">total ndofs: %d\n",sm->ndofs);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Use this to update ivar pos in each mat physics
    //++++++++++++++++++++++++++++++++++++++++++++++
    smat_physics_set_ivar_loc_array(sm->mat_physics_elem,sm->nmat_physics,&(sm->ivar_pos));

    //FORNOW JUST WORKS IN SERIAL
    #ifndef _MESSG
        sm->my_ndofs = sm->ndofs;
    #endif
    sm->ndofs_old = 0;
    sm->my_ndofs_old = 0;
    sm->sol       = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_old   = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_older = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->dirichlet_data = (double*) tl_alloc(sizeof(double), sm->ndofs);
    sm->bc_mask = (int*) tl_alloc(sizeof(int), sm->ndofs);

    sarray_init_dbl(sm->sol,sm->ndofs);
    sarray_init_dbl(sm->sol_old,sm->ndofs);
    sarray_init_dbl(sm->sol_older,sm->ndofs);
    sarray_init_dbl(sm->dirichlet_data,sm->ndofs);
    sarray_init_int(sm->bc_mask,sm->ndofs);


    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocating dependent variables
    //++++++++++++++++++++++++++++++++++++++++++++++
    smodel_super_build_dvars(sm);

    
    //free up anything local to routine
    ip = tl_free(sizeof(SIVAR_POSITION), grid->nnodes, ip);
    
}

