/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Intializes an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod (SMODEL_DESIGN *)  a pointer to an AdH design-level model
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_defaults(SMODEL_DESIGN *dm) {

    strcpy(dm->filename_design,"");
    strcpy(dm->filename_grid,"");
    strcpy(dm->filename_params,"");
    dm->nSuperModels = UNSET_INT;
    dm->superModel = NULL;
    dm->superCouple = NULL;
    dm->grid = NULL;
    dm->series_dt = NULL;
    dm->series_out= NULL;
    dm->params = NULL;
    dm->nUnique = UNSET_INT;
    dm->unique_id = NULL;
    dm->lin_sys_id = NULL;
    dm->lin_sys = NULL;

}
