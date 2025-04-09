/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This fil initializes and SIVAR_POS structure */
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
    ip->h = UNSET_INT;    strcpy(ip->name.h,"h");
    ip->u = UNSET_INT;    strcpy(ip->name.u,"u");
    ip->v = UNSET_INT;    strcpy(ip->name.v,"v");
    ip->w = UNSET_INT;    strcpy(ip->name.w,"w");
    ip->uda = UNSET_INT;  strcpy(ip->name.uda,"uda");
    ip->vda = UNSET_INT;  strcpy(ip->name.vda,"vda");
    ip->dpl = UNSET_INT;  strcpy(ip->name.dpl,"dpl");
    ip->prs = UNSET_INT;  strcpy(ip->name.prs,"prs");
    ip->heat = UNSET_INT; strcpy(ip->name.heat,"heat");
    ip->sal = UNSET_INT;  strcpy(ip->name.sal,"sal");
    sarray_init_value_int(ip->con,MAX_TRNS_VARS,UNSET_INT);
    for (int itrns=0; itrns<MAX_TRNS_VARS; itrns++) {
        char str[5]; sprintf(str, "%d", itrns);
        strcpy(ip->name.con[itrns],"con"); strcat(ip->name.con[itrns],str);
    }
}

void sivar_position_map(SIVAR_POSITION *ip, int *FLAG) {

    ip->n = 0;
    ip->ntrns = 0;

    if (FLAG[0] == 1) {ip->h    = ip->n; ip->n++;}
    if (FLAG[1] == 1) {ip->u    = ip->n; ip->n++;}
    if (FLAG[2] == 1) {ip->v    = ip->n; ip->n++;}
    if (FLAG[3] == 1) {ip->w    = ip->n; ip->n++;}
    if (FLAG[4] == 1) {ip->uda  = ip->n; ip->n++;}
    if (FLAG[5] == 1) {ip->vda  = ip->n; ip->n++;}
    if (FLAG[6] == 1) {ip->dpl  = ip->n; ip->n++;}
    if (FLAG[7] == 1) {ip->prs  = ip->n; ip->n++;}
    if (FLAG[8] == 1) {ip->heat = ip->n; ip->n++;}
    if (FLAG[9] == 1) {ip->sal  = ip->n; ip->n++;}
    for (int itrns=0; itrns<MAX_TRNS_VARS; itrns++) {
        if (FLAG[MAX_VARS + itrns] == 1) {
            ip->con[itrns] = ip->n; 
            ip->n++; 
            ip->ntrns++;
        }
    }

}

void sivar_position_printScreen(SIVAR_POSITION *ip) {

    printf("------ ip->n: %d\n",ip->n);
    printf("------ ip->ntrns: %d\n",ip->ntrns);
    printf("------ ip->%s: %d\n",ip->name.h,ip->h);
    printf("------ ip->%s: %d\n",ip->name.u,ip->u);
    printf("------ ip->%s: %d\n",ip->name.v,ip->v);
    printf("------ ip->%s: %d\n",ip->name.w,ip->w);
    printf("------ ip->%s: %d\n",ip->name.uda,ip->uda);
    printf("------ ip->%s: %d\n",ip->name.vda,ip->vda);
    printf("------ ip->%s: %d\n",ip->name.dpl,ip->dpl);
    printf("------ ip->%s: %d\n",ip->name.prs,ip->prs);
    printf("------ ip->%s: %d\n",ip->name.heat,ip->heat);
    printf("------ ip->%s: %d\n",ip->name.sal,ip->sal);
    for (int itrns=0; itrns<ip->ntrns; itrns++) {
        printf("------ ip->%s: %d\n",ip->name.con[itrns],ip->con[itrns]);
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
        if (ipNode->h != UNSET_INT) {
            iv[ip->h][i] = ndofs; ndofs++;
        }
        if (ipNode->u != UNSET_INT) {
            iv[ip->u][i] = ndofs; ndofs++;
        }
        if (ipNode->v != UNSET_INT) {
            iv[ip->v][i] = ndofs; ndofs++;
        }
        if (ipNode->w  != UNSET_INT) {
            iv[ip->w][i] = ndofs; ndofs++;
        }
        if (ipNode->uda != UNSET_INT) {
            iv[ip->uda][i]  = ndofs; ndofs++;
        }
        if (ipNode->vda != UNSET_INT) {
            iv[ip->vda][i]  = ndofs; ndofs++;
        }
        if (ipNode->dpl != UNSET_INT) {
            iv[ip->dpl][i]  = ndofs; ndofs++;
        }
        if (ipNode->prs != UNSET_INT) {
            iv[ip->prs][i]  = ndofs; ndofs++;
        }
        if (ipNode->heat != UNSET_INT) {
            iv[ip->heat][i] = ndofs; ndofs++;
        }
        if (ipNode->sal != UNSET_INT) {
            iv[ip->sal][i]  = ndofs; ndofs++;
        }
        for (int itrns=0; itrns<MAX_TRNS_VARS; itrns++) {
            if (ipNode->con[itrns] != UNSET_INT) {
                iv[ip->con[itrns]][i] = ndofs;
                ndofs++;
            }
        }
    }

    return ndofs;

}
