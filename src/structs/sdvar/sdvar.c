/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sdvar.c This file collects methods of the SCOVERAGE structure for coverages   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes a SDVAR instance
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dp (SDVAR_POSITION *)  a pointer to a SDVAR_POSITION structure

 * \note CJT - lists all possible AdH variables
 * \note CJT -  some of these variables are internally needed, others are turned on by user in paramter file, so it must be read first
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sdvar_position_init(SDVAR_POSITION *dp) {
    dp->n = 0;
    dp->var = (int *) tl_alloc(sizeof(int),adh_def.n_dvars);
    dp->print_flag = (bool *) tl_alloc(sizeof(bool),adh_def.n_dvars);
    
    for (int i=0; i<adh_def.n_dvars; i++) {
        dp->var[i] = UNSET_INT;
        dp->print_flag[i] = false;
    }
}

void sdvar_position_free(SDVAR_POSITION *dp) {
    dp->var = (int *) tl_free(sizeof(int),adh_def.n_dvars,dp->var);
    dp->print_flag = (bool *) tl_free(sizeof(bool),adh_def.n_dvars,dp->print_flag);
}

void sdvar_position_printScreen(SDVAR_POSITION *dp) {
    printf("SDVAR_POSITION ------------\n");
    printf("---- n: %d\n",dp->n);
    for (int i=0; i<adh_def.n_dvars; i++) {
        printf("---- Dependent Variable[%d]: %s - %s || position: %d || print_flag: %d\n",i,adh_def.dvar[i].name,adh_def.dvar[i].subname,dp->var[i],dp->print_flag[i]);
        
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and initializes a SDVAR dependent variable structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sdvar (SDVAR *)  a pointer to a SDVAR structure
 * @param[in] nnode (int)  number of active nodes
 * @param[in] n_dvar (int)  number of active dependent variables
 * @param[in] n_dvar_elem_dbl (int) the number of active elemental double variables
 * @param[in] n_dvar_elem_int (int) the number of active elemental integer variables
 * @param[in] nnode_dvar (int) the number of active nodes
 * @param[in] nelem_dvar (int) the number of active elements
 *
 * \note CJT -  some of these variables are internally needed, others are turned on by user in paramter file, so it must be read first
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sdvar_alloc_init(SDVAR *sdvar, int nnode, int n_dvar, int n_dvar_elem_dbl, int n_dvar_elem_int, int nnode_dvar, int nelem_dvar){
    int j;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++
    // initialize variable index positions and names
    // ++++++++++++++++++++++++++++++++++++++++++++++
    if (n_dvar > 0)          {sdvar_position_init(&sdvar->sdvar_pos_node);}
    if (n_dvar_elem_dbl > 0) {sdvar_position_init(&sdvar->sdvar_pos_elem_dbl);}
    if (n_dvar_elem_int > 0) {sdvar_position_init(&sdvar->sdvar_pos_elem_int);}
    
    // ++++++++++++++++++++++++++++++++++++++++++++++
    // allocate and init node-based variables
    // ++++++++++++++++++++++++++++++++++++++++++++++
    sdvar->n_dvar = n_dvar; //number of active dependent variables
    sdvar->nnode_dvar = nnode_dvar; //number of active nodes
    sdvar->nodal_dvar = NULL;
    if (sdvar->n_dvar > 0){
        sdvar->nodal_dvar = allocate_dptr_dbl(n_dvar,nnode_dvar);
        sarray_init_double_2d(sdvar->nodal_dvar, n_dvar, nnode_dvar);
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++
    // allocate and init element-based variables
    // ++++++++++++++++++++++++++++++++++++++++++++++
    sdvar->n_dvar_elem_dbl = n_dvar_elem_dbl;
    sdvar->nelem_dvar = nelem_dvar; //number of active elements
    sdvar->elem_dvar = NULL;
    if(sdvar->n_dvar_elem_dbl > 0){
        sdvar->elem_dvar = allocate_dptr_dbl(n_dvar_elem_dbl,nelem_dvar);
        sarray_init_double_2d(sdvar->elem_dvar,n_dvar_elem_dbl,nelem_dvar);
    }

    sdvar->n_dvar_elem_int = n_dvar_elem_int;
    sdvar->elem_flags = NULL;
    if(sdvar->n_dvar_elem_int > 0){
        sdvar->elem_flags = allocate_dptr_int(n_dvar_elem_int,nelem_dvar);
        sarray_init_int_2d(sdvar->elem_flags,n_dvar_elem_int,nelem_dvar);
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++
    // allocate and init maps
    // ++++++++++++++++++++++++++++++++++++++++++++++
    //maps, only need one if nnodes != nnode_dvar
    sdvar->dvar_node_map = NULL; //[nnode] array goes NodeID->index within nodal_dvar
    //if (nnode!=nnode_dvar){
        //create map here?, will also need the nodal physics to get this actually
    //}
    //same for dvar_active
    sdvar->dvar_active = NULL;
    //int *dvar_active; //n_dvar array that has the NodeID, convenient for printing out stuff
    //and for elem_map
    sdvar->dvar_elem_map = NULL;
    //int *dvar_elem_map; //nelem array takes elem # -> index within elem_var
}

