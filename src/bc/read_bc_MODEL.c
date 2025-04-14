/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file read_bc.c This file reads the MODEL card in a superfile input parameter file                                                              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void model_flag_switch(char *model, int *flags) {
    for (int i=0; i<adh_def.nmodels; i++) {
        if (strcmp(model,adh_def.model[i].subname) == 0) {flags[i] = 1;}
    }
}
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
 * \note MODEL [IMAT] [PHYSICS] ...
 * \note BC         [IMAT] [PHYSICS] [NB] [VEL] ...
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_bc_MODEL(SMODEL_SUPER *sm, FILE *fp) {

    int i,j,k,ivar;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token, str[MAXLINE] = "";
    int smat_id = UNSET_INT, imod = UNSET_INT, ipde = UNSET_INT;
    SIVAR_POSITION *ip = NULL;
    SMAT_PHYSICS *mat = NULL;
    SPDE *pde = NULL;
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading param file for MODEL and BC extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    
    //+++++++++++++++++++++++++++++++++++
    // first read to count materials
    //+++++++++++++++++++++++++++++++++++
    rewind(fp);
    sm->nmat_physics = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        // both MODEL & BC *CAN* create a new material
        if (strcmp(token, "MODEL") == 0 || strcmp(token, "BC") == 0) {
            smat_id = get_next_token_int(&token)-1;
            if (smat_id >= sm->nmat_physics) {sm->nmat_physics++;}
        }
    }
    
    //+++++++++++++++++++++++++++++++++++
    // allocated smat_physics array
    //+++++++++++++++++++++++++++++++++++
    smat_physics_alloc_init_array(&(sm->mat_physics_elem),sm->nmat_physics);
    
    //+++++++++++++++++++++++++++++++++++
    // read smat->id and smat->npdes
    //+++++++++++++++++++++++++++++++++++
    rewind(fp);
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        if (strcmp(token, "MODEL") == 0 || strcmp(token, "BC") == 0) {
            smat_id = get_next_token_int(&token)-1;
            mat = &(sm->mat_physics_elem[smat_id]);
            mat->id = smat_id;
            mat->npdes++;
        }
    }
    
    //+++++++++++++++++++++++++++++++++++
    // allocate smat pdes
    //+++++++++++++++++++++++++++++++++++
    for (i=0; i<sm->nmat_physics; i++) {
        mat = &(sm->mat_physics_elem[i]);
        printf("material[%d] || npdes: %d\n",i,mat->npdes);
        if (mat->npdes < 1) {
            sprintf(str, "ERROR: Material %d has %d physics residuals attached.\n",i,mat->npdes);
            tl_error(str);
        }
        mat->pde = (SPDE *) tl_alloc(sizeof(SPDE), mat->npdes);
        for (j=0; j<mat->npdes; j++) {
            mat->pde[j].imod = UNSET_INT;
            mat->pde[j].iseries = UNSET_INT;
            mat->pde[j].resid = NULL;
            mat->pde[j].init  = NULL;
        }
    }
    
    
    //+++++++++++++++++++++++++++++++++++
    // read SPDE data
    //+++++++++++++++++++++++++++++++++++
    rewind(fp);
    int ntrns_mod = 0, ntrns_bc = 0, series_id = UNSET_INT;
    int counts[sm->nmat_physics]; sarray_init_int(counts,sm->nmat_physics);
    char modType[MAXLINE], physType[MAXLINE], bcType[MAXLINE], varType[MAXLINE];
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue; // Read First Card
        strcpy(modType,token);
        
        smat_id = get_next_token_int(&token)-1; // Read Material ID
        mat = &(sm->mat_physics_elem[smat_id]);

        if (strcmp(modType, "MODEL") == 0) {
            get_next_token(&token); if (token == NULL) continue; // Read Physics Type
            model_flag_switch(token,sm->flags.model);
            pde = &mat->pde[counts[smat_id]];
            pde->resid = select_resid_func(token,&(pde->imod),ntrns_mod);
            pde->init = select_init_func(token,ntrns_mod);
            if (strcmp(token, "TRNS") == 0) {ntrns_mod++;}
            if (DEBUG) {printf("MODEL found: %s on material: %d || pde->imod: %d\n",token,smat_id,pde->imod);}
            counts[smat_id]++;
        } else if (strcmp(modType, "BC") == 0) {
            pde = &mat->pde[counts[smat_id]];
            get_next_token(&token); if (token == NULL) continue; // Read Physics Type
            strcpy(physType,token);
            get_next_token(&token); if (token == NULL) continue; // Read BC Type
            strcpy(bcType,token);
            get_next_token(&token); if (token == NULL) continue; // Read Variable Type
            strcpy(varType,token);
            pde->iseries = get_next_token_int(&token)-1; // Read Series ID
            if (strcmp(physType, "TRNS") == 0) {ntrns_bc++;}
            pde->resid = select_bresid_func(physType,bcType,varType,&(pde->imod),ntrns_bc);
            if (DEBUG) {printf("BC found: %s on material: %d || pde->imod: %d\n",token,smat_id,pde->imod);}
            counts[smat_id]++;
        }
    }
    rewind(fp);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Finished Reading
    // Use pde info to build up the physics materials
    // maps contained in mat->ivar_pos
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    smat_physics_set_ivar_pos(sm->mat_physics_elem, sm->nmat_physics);

}
