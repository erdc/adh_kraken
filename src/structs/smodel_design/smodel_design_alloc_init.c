/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL **)  a pointer to an AdH design-level model
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 * @param[in]  nSubModels                (int*) the total number of submodels in each supermodel
 * @param[in]  nFluxInterfaces      (int*) the total number of flux interfaces between supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_alloc_init(SMODEL_DESIGN **dmod, int nSuperModels) {
    int i;
    
    // allocate
    (*dmod) = (SMODEL_DESIGN *) tl_alloc(sizeof(SMODEL_DESIGN), 1);

    // initialize the design model and all its super and sub models
    (*dmod)->nSuperModels = nSuperModels;
    smodel_super_alloc_init(&((*dmod)->superModel),nSuperModels); //,nFluxInterfaces);

}
