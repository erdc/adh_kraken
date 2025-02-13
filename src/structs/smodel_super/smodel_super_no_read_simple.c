/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super_no_read_simple.c This file collects methods of the SUPER_MODEL structure */
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
void smodel_super_no_read_simple(SMODEL_SUPER *sm, double* dt_in, double* t_init, double* t_prev,
    double* t_final, int nmat_physics, char elemVarCode[4], int isSimple, SGRID *grid, SLIN_SYS *sys) {

    assert(nmat_physics == 1); // only 1 material for now

    int i,j;
    
    printf("Initializing smodel super without file read\n");
    // assign scalars
    sm->ts->dt = *dt_in;
    sm->ts->dt_old = *dt_in;
    sm->ts->dt_err = *dt_in;
    //sm->ts->dt_prev = *dt_in;
    sm->ts->t_init = *t_init;
    sm->ts->t_prev = *t_prev;
    sm->ts->t_final = *t_final;
    sm->nsubsteps = 1;
    sm->nseries = 0;
    sm->itrns = 0;
    sm->isSimple = isSimple;
    sm->grid = grid;
    sm->lin_sys = sys;

    sm->inc_nonlin = 1e-3;
    sm->tol_nonlin = 1e-5;

    // Set physics materials
    sm->nmat_physics = nmat_physics;
    if (sm->grid->nelems3d>0) {
        sm->elem3d_physics_mat = (int*) tl_alloc(sizeof(int), sm->grid->nelems3d);
        sarray_init_int(sm->elem3d_physics_mat, sm->grid->nelems3d); // only 1 material for now
    }
    if (sm->grid->nelems2d>0) {
        sm->elem2d_physics_mat = (int*) tl_alloc(sizeof(int), sm->grid->nelems2d);
        sarray_init_int(sm->elem2d_physics_mat, sm->grid->nelems2d); // only 1 material for now
    }
    if (sm->grid->nelems1d>0) {
        sm->elem1d_physics_mat = (int*) tl_alloc(sizeof(int), sm->grid->nelems1d);
        sarray_init_int(sm->elem1d_physics_mat, sm->grid->nelems1d); // only 1 material for now
    }
    char code[1][10]; strcpy(code,elemVarCode);
    smat_physics_alloc_init_array(&sm->mat_physics_elem, sm->nmat_physics,code);

    // nodal physics set-up
    sm->mat_physics_node = (SMAT_PHYSICS **) tl_alloc(sizeof(SMAT_PHYSICS *), sm->grid->nnodes);
    for(i=0;i<sm->grid->nnodes;i++){
        sm->mat_physics_node[i] = &(sm->mat_physics_elem[0]); // only 1 material for now
    }


    int FLAGS[MAX_TRNS_VARS + MAX_VARS]; sarray_init_int(FLAGS,MAX_TRNS_VARS + MAX_VARS);
    sivar_position_init(&(sm->ivar_pos));
    smat_physics_position_flag(sm->mat_physics_node,grid->nnodes,FLAGS); 
    sivar_position_map(&(sm->ivar_pos),FLAGS);    

    SIVAR_POSITION ip[grid->nnodes];
    for (i=0; i<grid->nnodes; i++) {ip[i] = sm->mat_physics_node[i]->ivars;}
    sm->ndofs = sivar_position_build_dof_map(&sm->ivar_pos,grid->nnodes,ip,&sm->ivars);

    sm->forward_step = FE_NEWTON;


}
