/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_printScreen.c This file prints an  SMAT_PHYSICS structures          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Flags independent variables on all nodes to get all ivariables used on superModel
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] mat (SMAT_ELEM *)  pointer to a transport material
 * @param[in] id  (int) the SMAT ID
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_printScreen(SMAT_PHYSICS *m) {
    int itrns;

    printf("----------------------------------------\n");
    printf("SMAT #%d || code: %s\n",m->id+1,m->code);
    printf("----------------------------------------\n");
    sivar_position_printScreen(&m->ivar_pos);
    printf("------ SW_FLOW: %d \n",m->SW_FLOW);
    printf("------ SW1_FLOW: %d \n",m->SW1_FLOW);
    printf("------ SW2_FLOW: %d \n",m->SW2_FLOW);
    printf("------ SW3_FLOW: %d \n",m->SW3_FLOW);
    printf("------ NS_FLOW: %d \n",m->NS_FLOW);
    printf("------ NS3_FLOW: %d \n",m->NS3_FLOW);
    printf("------ NS3_SPLIT: %d \n",m->NS3_SPLIT);
    printf("------ DW_FLOW: %d \n",m->DW_FLOW);
    printf("------ WVEL_SPLIT: %d \n",m->WVEL_SPLIT);
    printf("------ PRESSURE: %d \n",m->PRESSURE);
    printf("------ GW_FLOW: %d \n",m->GW_FLOW);
    printf("------ VORTICITY: %d \n",m->VORTICITY);
    printf("------ SEDIMENT: %d \n",m->SEDIMENT);
    printf("------ SEDLIB: %d \n",m->SEDLIB);
    printf("------ ICM: %d \n",m->ICM);
    printf("------ NSM: %d \n",m->NSM);
    printf("------ WAVE: %d \n",m->WAVE);
    printf("------ WIND: %d \n",m->WIND);

    printf("------ ntrns: %d\n",m->ntrns);
    for (itrns=0; itrns<m->ntrns; itrns++) {
        printf("------ TRANSPORT[%d]: %d\n",itrns,m->TRANSPORT[itrns]);
    }

    printf("------ nsubmodels: %d\n",m->nSubmodels);
    for (int i=0; i<m->nSubmodels; i++) {
        smodel_printScreen(&m->model[i]);
    }
    //printf("----------------------------------------\n");

}
