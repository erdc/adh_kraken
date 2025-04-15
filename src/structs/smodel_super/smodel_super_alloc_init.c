/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an array of AdH Super Models
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod           (SUPER_MODEL **)  a double pointer to an array of AdH supermodels
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 *
 * \note inherets an SGRID and SDT from designer mode
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_alloc_init(SMODEL_SUPER **superModel, int nSuperModels) {
    
    // assertions
    assert(nSuperModels > 0);  // must have a least one superModel defined to run
    
    // allocate all the supermodels
    (*superModel) = (SMODEL_SUPER *) tl_alloc(sizeof(SMODEL_SUPER), nSuperModels);

    // initialize
    SMODEL_SUPER *sm;
    for (int isup=0; isup<nSuperModels; isup++) {
        sm = &((*superModel)[isup]); // alias
        sflags_init(&(sm->flags));
        strcpy(sm->filebase,"");

#ifdef _MESSG
        sm->myid = messg_comm_rank(sm->grid->smpi->ADH_COMM);
        sm->npes = messg_comm_size(sm->grid->smpi->ADH_COMM);
#else
        sm->myid = 0;
        sm->npes = 1;
#endif
        sm->itrns = UNSET_INT;
        sm->isSimple = UNSET_INT;
        sm->id = isup;
        sm->type = NULL;
        sm->grid = NULL;
        sm->ts = NULL;
        sm->nsubsteps = 1;
        sm->meshcode = 0;
        //Eventually need to hook into front end
        //but only FE_NEWTON time step available right now
        sm->forward_step = FE_NEWTON;
        sm->sw = NULL;
        // sm->con = NULL;
        // sm->sed = NULL;
        // sm->gw = NULL;

        //++++++++++++++++++++++++++++++++++++++++
        // linear solver defaults
        //++++++++++++++++++++++++++++++++++++++++
        sm->lin_sys = NULL;
        sm->perturbation  = sqrt(SMALL);
        sm->inc_nonlin = 1e-8;
        sm->tol_nonlin = 1e-10;
        sm->max_nonlin_linesearch_cuts = 5;
        sm->max_nonlin_it = 10;
        sm->it_count_nonlin = 0;
        sm->force_nonlin_it = NO;
        sm->LINEAR_PROBLEM = NO;
        sm->it_count_nonlin_failed = 0;
        sm->nalloc_inc = 1;
        sm->nonlinear_it_total = 0;
        //++++++++++++++++++++++++++++++++++++++++

        sm->my_ndofs = UNSET_INT;
        sm->my_ndofs_old = UNSET_INT;
        sm->ndofs = UNSET_INT;
        sm->ndofs_old = UNSET_INT;
        sm->macro_ndofs = UNSET_INT;
        sm->macro_ndofs_old = UNSET_INT;
        sm->sol = NULL;
        sm->sol_old = NULL;
        sm->sol_older = NULL;
        sm->elem1d_physics_mat = NULL;
        sm->elem2d_physics_mat = NULL;
        sm->elem3d_physics_mat = NULL;
        sm->nmat_physics = 0;
        sm->mat_physics_elem = NULL;
        sm->mat_physics_node = NULL;
        sm->ivars = NULL;
        sm->bc_mask = NULL;
        sm->dirichlet_data = NULL;
        sm->gravity = 9.8;
        sm->density = 1000.;
        sm->o_flag = 0;
        sm->initial_grid_mass = 0.;
        sm->grid_mass_error = 0.;
        sm->nseries = 0;
        sm->series_head = NULL;
        sm->series_curr = NULL;
        sm->series_wind_head = NULL;
        sm->series_wind_curr = NULL;
        sm->series_wave_head = NULL;
        sm->series_wave_curr = NULL;
        sm->series_gw_psk_head = NULL;
        sm->series_gw_psk_curr = NULL;

    }   
}

