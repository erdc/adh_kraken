/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_printScreen.c This file prints an SMODEL struct to the screen             */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     prints an SMODEL struct to the screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] model (SMODEL **)  the struct double pointer
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_printScreen(SMODEL *m) {
    printf("-------- model physics ID: %d\n",m->physics); // a model ID for elemental residual functions, etc.
    printf("-------- model physics init ID: %d\n",m->physics_init);   // a model ID for elemental body initialization
    printf("-------- model nvariables: %d\n",m->nvar);           // the # of independent variable on this model
    for (int i=0; i<m->nvar; i++) {
        printf("-------- model indepedent variable %d position: %d\n",i,m->physics_vars[i]);
    }
}