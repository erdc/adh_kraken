#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_free_array(SMODEL_SUPER *sm, int nSuper) {
    int i;
    
    for (i=0; i<nSuper; i++) {smodel_super_free(&(sm[i]));}
    sm = (SMODEL_SUPER *) tl_free(sizeof(SMODEL_SUPER), nSuper, sm);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_free(SMODEL_SUPER *sm) {


    int ntransport=0,nlayer=0,nsed=0; // CJT -- send this later

    if (sm->ts != NULL) {
        sm->ts = (SDT *) tl_free(sizeof(SDT), 1, sm->ts);
    }
        //if (sm->con != NULL) {scon_free(sm->con);}
        //if (sm->sed != NULL) {ssed_free(sm->sed);}
        //if (sm->gw != NULL) {sgw_free(sm->gw);}

    if (sm->sol != NULL) {
        sm->sol = (double *) tl_free(sizeof(double), sm->ndofs, sm->sol);
    }
    if (sm->sol_old != NULL) {
        sm->sol_old = (double *) tl_free(sizeof(double), sm->ndofs, sm->sol_old);
    }
    if (sm->sol_older != NULL) {
        sm->sol_older = (double *) tl_free(sizeof(double), sm->ndofs, sm->sol_older);
    }

    if (sm->elem1d_physics_mat != NULL) {
        sm->elem1d_physics_mat = (int *) tl_free(sizeof(int), sm->grid->nelems1d, sm->elem1d_physics_mat);
    }
    if (sm->elem2d_physics_mat != NULL) {
        sm->elem2d_physics_mat = (int *) tl_free(sizeof(int), sm->grid->nelems2d, sm->elem2d_physics_mat);
    }
    if (sm->elem3d_physics_mat != NULL) {
        sm->elem2d_physics_mat = (int *) tl_free(sizeof(int), sm->grid->nelems3d, sm->elem3d_physics_mat);
    }

    if (sm->mat_physics_elem != NULL) {smat_physics_free_array(sm->mat_physics_elem,1);}
    if (sm->mat_physics_node != NULL) {
        sm->mat_physics_node = (SMAT_PHYSICS **) tl_free(sizeof(SMAT_PHYSICS *), sm->grid->nnodes, sm->mat_physics_node);
    }

    if (sm->ivars != NULL) {
        for (int iv=0; iv<sm->ivar_pos.n; iv++) {
            sm->ivars[iv] = (int *) tl_free(sizeof(int), sm->grid->nnodes, sm->ivars[iv]);
        }
        sm->ivars = (int **) tl_free(sizeof(int *), sm->ivar_pos.n, sm->ivars);
    }

    if (sm->series_wind_head != NULL) {sseries_free(sm->series_wind_head);}
    if (sm->series_wind_curr != NULL) {sseries_free(sm->series_wind_curr);}

    if (sm->series_wave_head != NULL) {sseries_free(sm->series_wave_head);}
    if (sm->series_wave_curr != NULL) {sseries_free(sm->series_wave_curr);}


    if (sm->series_gw_psk_head != NULL) {sseries_free(sm->series_gw_psk_head);}
    if (sm->series_gw_psk_curr != NULL) {sseries_free(sm->series_gw_psk_curr);}

        //if (sm->str_values != NULL) {sstr_value_free(sm->str_values,ntransport,nlayer,nsed);}

    if(sm->bc_mask!=NULL){
        sm->bc_mask = (int *) tl_free(sizeof(int), (sm->ndofs), sm->bc_mask);
    }
    if(sm->dirichlet_data!=NULL){
        sm->dirichlet_data = (double *) tl_free(sizeof(double), (sm->ndofs), sm->dirichlet_data);
    }

}