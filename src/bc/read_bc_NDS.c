/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_NDS.c  This file collects methods to read node strings                                                                                 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel parameter file for node strings, allocates and fills in the double array
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] fp (FILE *) a superModel parameter file
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_NDS(SMODEL_SUPER *sm, FILE *fp) {
    
    int i;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token, str[MAXLINE];
    int nodeID = UNSET_INT, lnodeID = UNSET_INT, nnodes_string = 0;
    int *iarray = NULL;
    
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        if (strcmp(token, "NDS") == 0) {

            // check type of node string (can be boundary condition, print, etc.)
            get_next_token(&token);
            if (strcmp(token, "BC") == 0) {
                iarray = sm->node_strings_bc;
                iarray = (int *) tl_alloc(sizeof(int), sm->grid->nnodes);
                sarray_init_value_int(iarray,sm->grid->nnodes,UNSET_INT);
            }
            
            // get number of nodes in this string
            nnodes_string = get_next_token_int(&token);
            
            // read nodes
            for (i=0; i<nnodes_string; i++) {
                read = getline(&line, &len, fp);
                get_token(line,&token);
                sscanf(token, "%d", &nodeID);
                nodeID--;
                lnodeID = is_my_node(nodeID, sm->grid); // go from global to local ID
                if (lnodeID < 0){break;}                // global node not on this processor
                iarray[lnodeID]  = get_next_token_int(&token); // get material
            }
        }
    }
    rewind(fp);
}
