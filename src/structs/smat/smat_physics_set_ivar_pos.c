/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_set_ivar_pos.c The file is initializing an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     sets the ivar_pos within the smat_physics objects
 *             given the pde info is filled out
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
void smat_physics_set_ivar_pos(SMAT_PHYSICS *mat_physics_elem, int nmat_physics){
	int varID = UNSET_INT;
	SMAT_PHYSICS *mat = NULL;
	SIVAR_POSITION *ip = NULL;
	int ctr, imod;
    int FLAGS[adh_def.n_ivars];


	//Loop over each mat physics
    //We will always order
    //variables in same order as 
    //they are in adh_def
    for (int i=0; i<nmat_physics; i++) {
        sarray_init_value_int(FLAGS, adh_def.n_ivars, UNSET_INT);
        mat = &(mat_physics_elem[i]);
        ip = &mat->ivar_pos;
        // Loop over each pde within the mat
        for (int ipde=0; ipde<mat->npdes; ipde++) {
            imod = mat->pde[ipde].imod;  // alias, this is model # in adh_def
            for (int ivar=0; ivar<adh_def.model[imod].nivars; ivar++) {
                varID = adh_def.model[imod].var[ivar];
                // only add if variable hasn't been used in another subModel
                // if varID is repeated on a material this should result in 
                // error (Can't solve for same variable twice on one element)
                if (FLAGS[varID] != UNSET_INT) {
                	//error message
                	tl_error("Can only assign one pde per solution variable");
                }else{
                    FLAGS[varID] = 1;
                }
            }
        }
        //reallocate var_code so we dont carry around extra stuff\
        //retain first ip->n entries
        sivar_position_map(ip, FLAGS);
    

        // This will be set later, needs info for entire super model to set
        mat->ivar_loc = (int *) tl_alloc(sizeof(int), ip->n);
        sarray_init_value_int(mat->ivar_loc, mat->ivar_pos.n, UNSET_INT);
        //will be set in smat_physics_set_ivar_loc
        
    }
}