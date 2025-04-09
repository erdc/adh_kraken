/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_alloc_init.c This file collects methods to allocate an SMODEL structure   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array SMODEL structures
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mod        (SMODEL **) the struct double pointer
 * @param[in]    phys       (int)   a model ID for elemental residual functions, etc.
 * @param[in]    phys_init  (int)   a model ID for elemental body initialization
 * @param[in]    nmods      (int)   the number of SMODELs to allocate
 * @param[in]    nvars      (int *) the number of variables in each submodel
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_alloc_init_array(SMODEL **mod, int *phys, int *physInit, int nmods, int *nvars) { 
    assert(nmods>0);
    (*mod) = (SMODEL *) tl_alloc(sizeof(SMODEL), nmods);
    for (int i=0; i<nmods; i++) {
    	smodel_alloc_init( &(*mod)[i], phys[i], physInit[i], nvars[i]);
    }
    return;
}   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates an SMODEL structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_alloc_init(SMODEL *m, int phys, int physInit, int nvar) {
    assert(nvar>0);
    m->physics = phys;
    m->physics_init = physInit;
    m->nvar = nvar;
    (m->physics_vars) = (int*) tl_alloc(sizeof(int), nvar);
    for(int i=0; i<nvar; i++){
        m->physics_vars[i] = 0;
    }
    return;
}   