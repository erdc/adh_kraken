/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints an AdH Designer Model to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL *)  a pointer to an AdH design-level model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_printScreen(SMODEL_DESIGN *dmod) {
    int iSuperMod = 0;

    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("nSuperModels: %d\n",dmod->nSuperModels);
    for (iSuperMod=0; iSuperMod<dmod->nSuperModels; iSuperMod++) {
        printf("--------------------------------------------\n");
        printf("SuperModel #%d \n",iSuperMod);
        printf("--------------------------------------------\n");
        smodel_super_printScreen(&(dmod->superModel[iSuperMod]));
    }
    //sgrid_printScreen(dmod->grid);
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
}
