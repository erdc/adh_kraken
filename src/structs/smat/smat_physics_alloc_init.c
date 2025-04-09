/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_alloc_init The file is for allocating an initializing an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 * @param[in]    nmat        (int) number of materials to allocate
 * @param[in]    codes       (char [][10]) an array of strings that encode the materials physics
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat) {
    assert(nmat >= 0);
    if (nmat == 0){
        mat_physics = NULL;
        return;
    } else {
        (*mat_physics) = (SMAT_PHYSICS *) tl_alloc(sizeof(SMAT_PHYSICS), nmat);
        SMAT_PHYSICS *mat;  // alias
        for (int imat=0; imat<nmat; imat++) {
            mat = &(*mat_physics)[imat];
            mat->id = imat;
            smat_physics_alloc_init(mat);
        }
        return;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 * @param[in]    codes       (char [][10]) an array of strings that encode the materials physics
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_alloc_init(SMAT_PHYSICS *mat) {

    sivar_position_init(&mat->ivar_pos);
    mat->n = 0;
    mat->ntrns = 0;
    mat->ivar_loc = NULL;
    mat->npdes = 0;
    mat->pde = NULL;

    //Remainder of models will be filled in in smat_physics_update_array.c since it requires **ivar map
}

