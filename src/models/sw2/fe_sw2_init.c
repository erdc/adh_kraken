/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_init.c This file is the init routine for sw2
 *         which will be called right before each Newton solve       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Performs any preprocessing for dependent variables in SW2 equations.
 *  For now this mainly is updating the wet/dry flags and the dirichlet conditions
 *  for dry elements
 *  \author    Mark Loveland, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] sm (SMODEL_SUPER *) a super model
 *
 * \note CJT\:: Label wetting and drying AND neighboring elements as wet/dry elements
 * \note CJT\:: flag = 0 :: fully wet || flag = 1 :: some nodes are dry || flag = 2 :: all nodes dry
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//needs to be called inside an init routine
int fe_sw2_init(SMODEL_SUPER *sm) {
    fe_sw2_wdflag_legacy(sm);
    // update Dirichlet condition for nodes with depth less than zero
    //how are mappings going to look??
    //this should work for now
    int i,temp,tempu,tempv;
    int h_loc = sm->ivar_pos.var[adh_def._H];
    int u_loc = sm->ivar_pos.var[adh_def._UDA];
    int v_loc = sm->ivar_pos.var[adh_def._VDA];
    //only want to loop over active nodes here not all nodes!!!! Will lead to bug
    for (i=0; i<sm->grid->nnodes; i++) {
        //will only work for CG, will need to generalize
        temp = sm->ivars[h_loc][i];
        
        if (sm->sol_old[temp] <= 0.) {
            //get the two corresponding u and v entries
            tempu = sm->ivars[u_loc][i];
            tempv = sm->ivars[v_loc][i];
            sm->bc_mask[tempu] = YES;
            sm->bc_mask[tempv] = YES;
        }
    }

    
    return 0;
}
