/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL *)  a pointer to an AdH design-level model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_free(SMODEL_DESIGN *dm) {

    if (dm->params != NULL) {scoverage_free(dm->params);}
    if (dm->series_dt != NULL) {sseries_free(dm->series_dt);}
    if (dm->series_out != NULL) {sseries_free(dm->series_out);}

    if (dm->lin_sys_id != NULL) {
        dm->lin_sys_id = (int *) tl_free(sizeof(int), dm->nlin_sys, dm->lin_sys_id);
    }
    //slin_sys_free(dm->lin_sys);

    if (dm->superCouple != NULL) {
        dm->superCouple = (SMODEL_COUPLE *) tl_free(sizeof(SMODEL_COUPLE), dm->nSuperModels-1, dm->superCouple);
    }
    
    if (dm->nSuperModels != UNSET_INT){
        smodel_super_free_array(dm->superModel,dm->nSuperModels);
    }

    sgrid_free(dm->grid);
    //nFluxInterfaces = (int *) tl_free(sizeof(int), sm_p->nSuperModels, nFluxInterfaces);
    //dm = (SMODEL_DESIGN *) tl_free(sizeof(SMODEL_DESIGN), 1, dm);
}
