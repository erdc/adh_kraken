/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sivar_position.c This fil initializes and SIVAR_POS structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Intializes the SIVAR_POS structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] ip (SIVAR_POS *)  a pointer to a SIVAR POS structure
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sivar_position_init(SIVAR_POSITION *ip) {
    
    ip->n = 0;
    ip->ntrns = 0;
    sarray_init_value_int(ip->var,N_IVARS_TOTAL,UNSET_INT);
}

void sivar_position_map(SIVAR_POSITION *ip, int *FLAG) {

    ip->n = 0;
    ip->ntrns = 0;
    for (int i=0; i<N_IVARS_TOTAL; i++) {
        if (FLAG[i] == 1) {
            ip->var[i] = ip->n;
            ip->n++;
            if (i > N_IVARS - 1) {ip->ntrns++;}
        }
    }
}

void sivar_position_printScreen(SIVAR_POSITION *ip) {

    printf("------ ip->n: %d\n",ip->n);
    printf("------ ip->ntrns: %d\n",ip->ntrns);
    for (int i=0; i<N_IVARS_TOTAL; i++) {
        printf("-------- ip->%s: %d\n",IVAR_NAME[i],ip->var[i]);
    }
    
}

int sivar_position_build_dof_map(SIVAR_POSITION *ip, int nnodes, SIVAR_POSITION *ipNode, int ***ivars) {

    int i;

    (*ivars) = (int **) tl_alloc(sizeof(int *), ip->n);
    int **iv = *(ivars);

    for (i=0; i<ip->n; i++) {
        iv[i] = (int *) tl_alloc(sizeof(int), nnodes);
        sarray_init_value_int(iv[i],nnodes,UNSET_INT);
    }
    
    int ndofs = 0;
    for (i=0; i<nnodes; i++) {
        for (int j=0; j<N_IVARS_TOTAL; j++) {
            if (ipNode->var[j] != UNSET_INT) {
                iv[ip->var[j]][i] = ndofs; ndofs++;
            }
        }
    }

    return ndofs;

}
