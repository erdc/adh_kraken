/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file read_bc.c This file reads the MODEL card in a superfile input parameter file                                                              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supersmels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_bc_MODEL(SMODEL_SUPER *sm, FILE *fp, char codes[][10]) {

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token;
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading param file for MODEL extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    
    sm->nmat_physics = 0;
    int smat_id = UNSET_INT, itrns = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token);
        if (token == NULL) continue;
        if (strcmp(token, "MODEL") == 0) {
            smat_id = get_next_token_int(&token)-1;
            if (smat_id >= sm->nmat_physics) {sm->nmat_physics++;}

            get_next_token(&token); if (token == NULL) continue;
            if      (strcmp(token,"SW1")  == 0) {codes[smat_id][0] = '1';} 
            else if (strcmp(token,"SW2")  == 0) {codes[smat_id][0] = '2';} 
            else if (strcmp(token,"SW3")  == 0) {codes[smat_id][0] = '3';} 
            else if (strcmp(token,"NS")   == 0) {codes[smat_id][0] = '4';} 
            else if (strcmp(token,"NSP")  == 0) {codes[smat_id][0] = '5';} 
            else if (strcmp(token,"DW")   == 0) {codes[smat_id][0] = '6';}  
            else if (strcmp(token,"WVEL") == 0) {codes[smat_id][0] = '7';} 
            else if (strcmp(token,"PRS")  == 0) {codes[smat_id][0] = '8';} 
            else if (strcmp(token,"GW")   == 0) {codes[smat_id][1] = '1';} 
            else if (strcmp(token,"TRNS") == 0) {
                itrns++;
                assert(itrns < 16);
                if (itrns == 1) codes[smat_id][3] = '1';
                if (itrns == 2) codes[smat_id][3] = '2';
                if (itrns == 3) codes[smat_id][3] = '3';
                if (itrns == 4) codes[smat_id][3] = '4';
                if (itrns == 5) codes[smat_id][3] = '5';
                if (itrns == 6) codes[smat_id][3] = '6';
                if (itrns == 7) codes[smat_id][3] = '7';
                if (itrns == 8) codes[smat_id][3] = '8';
                if (itrns == 9) codes[smat_id][3] = '9';
                if (itrns == 10) {codes[smat_id][2] = '1'; codes[smat_id][3] = '0';}
                if (itrns == 11) {codes[smat_id][2] = '1'; codes[smat_id][3] = '1';}
                if (itrns == 12) {codes[smat_id][2] = '1'; codes[smat_id][3] = '2';}
                if (itrns == 13) {codes[smat_id][2] = '1'; codes[smat_id][3] = '3';}
                if (itrns == 14) {codes[smat_id][2] = '1'; codes[smat_id][3] = '4';}
                if (itrns == 15) {codes[smat_id][2] = '1'; codes[smat_id][3] = '5';}
            } 
            else {tl_error("Model within SuperModel not recognized.\n");}

            //printf("code[smat_id = %d]: %d\n",smat_id,codes[smat_id]);
        }
    }
}
