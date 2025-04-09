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
void smodel_super_no_read(SMODEL_SUPER *sm, char **codes, int nmat_physics, int *mat_ids) {
    int i;
    SGRID *grid = sm->grid; // alias
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Set Element Physics Materials with codes input
    //++++++++++++++++++++++++++++++++++++++++++++++
    sm->nmat_physics = nmat_physics;

    //smat_physics_alloc_init_array(&(sm->mat_physics_elem),sm->nmat_physics,codes);
    smat_physics_alloc_init_array(&(sm->mat_physics_elem),sm->nmat_physics);
    if (DEBUG) {
        for (i=0; i<sm->nmat_physics; i++) {smat_physics_printScreen(&sm->mat_physics_elem[i]);}
    }
    
    
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
    int node_nvars[grid->nnodes]; sarray_init_int(node_nvars,grid->nnodes);
    
    int ie,imat,nd;
    for (ie=0; ie<grid->nelems1d; ie++) {
        if (grid->elem1d[ie].bflag != BODY) continue;
        imat = sm->elem1d_physics_mat[ie];
        for (i=0; i<grid->elem1d[ie].nnodes; i++) {
            nd = grid->elem1d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem1d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems2d; ie++) {
        if (grid->elem2d[ie].bflag != BODY) continue;
        imat = sm->elem2d_physics_mat[ie];
        for (i=0; i<grid->elem2d[ie].nnodes; i++) {
            nd = grid->elem2d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem2d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems3d; ie++) {
        imat = sm->elem3d_physics_mat[ie];
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            nd = grid->elem3d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem3d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    
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
    printf(">Calling smat_physics_update_array\n");
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Use this to update ivar pos in mat physics
    // also sets proper models
    //++++++++++++++++++++++++++++++++++++++++++++++
    smat_physics_update_array(sm->mat_physics_elem,sm->nmat_physics,&sm->ivar_pos);

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Use this info to fill in resid pointers
    // if we want to go object-oriented way
    //++++++++++++++++++++++++++++++++++++++++++++++
    //smodel_super_set_resid_init(sm, sm->mat_physics_elem, sm->nmat_physics);

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate the map from the independent variables 
    // to their position in the solution vector
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">allocating solution variables\n");
    }
    //this wont scale, need to dynammically allocate
    //SIVAR_POSITION ip[grid->nnodes];
    SIVAR_POSITION *ip;
    ip = (SIVAR_POSITION *) tl_alloc(sizeof(SIVAR_POSITION), grid->nnodes);
    
    for (i=0; i<grid->nnodes; i++) {ip[i] = sm->mat_physics_node[i]->ivar_pos;}
    sm->ndofs = sivar_position_build_dof_map(&sm->ivar_pos,grid->nnodes,ip,&sm->ivars);
    //FORNOW JUST WORKS IN SERIAL
#ifndef _MESSG
    sm->my_ndofs = sm->ndofs;
    sm->macro_ndofs = sm->ndofs;
#endif
    sm->ndofs_old = 0;
    sm->my_ndofs_old = 0;
    sm->macro_ndofs_old = 0;

    if (DEBUG) {
        printf(">total ndofs: %d\n",sm->ndofs);
    }
    sm->sol       = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_old   = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_older = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->dirichlet_data = (double*) tl_alloc(sizeof(double), sm->ndofs);
    sm->bc_mask = (int*) tl_alloc(sizeof(int), sm->ndofs);
    
    sm->tol_nonlin = 1e-7;
    sm->inc_nonlin = 1e-5;
    sm->max_nonlin_linesearch_cuts = 5;
    sm->it_count_nonlin_failed = 0;
    sm->max_nonlin_it = 25;
    sm->LINEAR_PROBLEM = NO;
    sm->force_nonlin_it = NO;
    sm->force_nonlin_it = NO;
    sm->nonlinear_it_total = 0;


    sarray_init_dbl(sm->sol,sm->ndofs);
    sarray_init_dbl(sm->sol_old,sm->ndofs);
    sarray_init_dbl(sm->sol_older,sm->ndofs);
    sarray_init_dbl(sm->dirichlet_data,sm->ndofs);
    sarray_init_int(sm->bc_mask,sm->ndofs);

    //free up anything local to routine
    ip = tl_free(sizeof(SIVAR_POSITION), grid->nnodes, ip);
    
}

