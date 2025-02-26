/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super_read.c This file reads user input data for a SUPER_MODEL structure  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads in SuperModel data
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm  (SUPER_MODEL *)  pointer to an AdH superModel
 * @param[in]  g      (SGRID *) a designer level grid the superModel resides on
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_read(SMODEL_SUPER *sm) {
    
    int i;
    SGRID *grid = sm->grid; // alias
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token,mode[MAXLINE],ext[MAXLINE];
    SFILE file;
    FILE *fp = NULL;

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Open SuperModel File
    //++++++++++++++++++++++++++++++++++++++++++++++
    strcpy(mode,"r"); strcpy(ext,".params");
    sfile_open(&file,sm->filebase,NULL,NULL,ext,mode,TRUE);
    fp = file.fp;
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read Model Parameters
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">reading param file for series extraction\n");
    }
    int type = UNSET_INT;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        if (strcmp(token, "SERIES") == 0) {
            get_next_token(&token);
            if (strcmp(token, "WIND") == 0) {
                sseries_read_allocate(sm->series_wind_head,sm->series_wind_curr,WIND_SERIES,&token,UNSET_INT);
            } else if (strcmp(token, "WAVE") == 0) {
                sseries_read_allocate(sm->series_wave_head,sm->series_wave_curr,WAVE_SERIES,&token,UNSET_INT);
            } else if (strcmp(token, "GWCONSTITUENT") == 0) {
                //sseries_read_allocate(sm->series_gw_psk_head,sm->series_gw_psk_curr,GWCONSTITUENT_SERIES,&token,UNSET_INT);
            } else {
                tl_error("Series not recognized.\n");
            }
        }
    }
    rewind(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read Element Physics Materials
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">reading bc file for model extraction\n");
    }
    // this needs to get moved to strings
    char codes[50][10]; 
    for (i=0; i<50; i++) {strcpy(codes[i],"0000");}
    read_bc_MODEL(sm,fp,codes);
    if (DEBUG) {
        printf("Element Material Info:\n");
        printf("-- Number of physics materials found: %d\n",sm->nmat_physics);
        for (i=0; i<sm->nmat_physics; i++) {
            printf("-- physics code[%d]: %s\n",i,codes[i]);
        }
    }
    smat_physics_alloc_init_array(&(sm->mat_physics_elem),sm->nmat_physics,codes);
    rewind(fp);
    if (DEBUG) {
        for (i=0; i<sm->nmat_physics; i++) {smat_physics_printScreen(&sm->mat_physics_elem[i]);}
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    //tl_check_all_pickets(__FILE__,__LINE__);//exit(-1);
    if (DEBUG) {
        printf(">reading bc file for string extraction\n");
    }
    while ( (read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token);
        if (token == NULL) continue;
        if (strcmp(token, "NDS") == 0 ||
            strcmp(token, "EGS") == 0 ||
            strcmp(token, "MDS") == 0 ||
            strcmp(token, "MTS") == 0) {
            //read_bc_STRINGS(sm,&token);
        }
    }
    rewind(fp);
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">reading rest of bc file\n");
    }
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        //  if      (strcmp(token, "OP")  == 0) {read_bc_OP(sm,&token);}
        //  else if (strcmp(token, "IP")  == 0) {read_bc_IP(sm,&token);}
        //  else if (strcmp(token, "PC")  == 0) {read_bc_PC(sm,&token);}
        //  else if (strcmp(token, "TC")  == 0) {read_bc_TC(sm,&token);}
        //  else if (strcmp(token, "OFF") == 0) {read_bc_OFF(sm,&token);}
        //  else if (strcmp(token, "MP")  == 0) {read_bc_MP(sm,&token);}
        //  else if (strcmp(token, "DB")  == 0) {read_bc_DB(sm,&token);}
        //  else if (strcmp(token, "NB")  == 0) {read_bc_NB(sm,&token);}
        //  else if (strcmp(token, "CN")  == 0) {read_bc_CN(sm,&token);}
        // //else if (strcmp(token, "SEDLIB") == 0) {read_bc_SEDLIB(sm,line);}
        //  else if (strcmp(token, "OP") == 0) {read_bc_OP(sm,&token);}
    }
    rewind(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    fclose(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Elemental Physics Coverage File
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">reading physics coverage file extraction\n");
    }
    strcpy(mode,"r"); strcpy(ext,".physics");
    sfile_open(&file,sm->filebase,NULL,NULL,ext,mode,TRUE);
    fp = file.fp;
    if (sm->grid->nelems3d > 0) {sm->elem3d_physics_mat = (int *) tl_alloc(sizeof(int),sm->grid->nelems3d);}
    if (sm->grid->nelems2d > 0) {sm->elem2d_physics_mat = (int *) tl_alloc(sizeof(int),sm->grid->nelems2d);}
    if (sm->grid->nelems1d > 0) {sm->elem1d_physics_mat = (int *) tl_alloc(sizeof(int),sm->grid->nelems1d);}
    int n1d=0,n2d=0,n3d=0,nmat1d[30],nmat2d[30],nmat3d[30];
    sarray_init_int(nmat1d,30); 
    sarray_init_int(nmat2d,30); 
    sarray_init_int(nmat3d,30);
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token);
        if (token == NULL) continue;
        if ( (strcmp(token,"TET") == 0) || (strcmp(token,"PRISM") == 0) ) {
            sm->elem3d_physics_mat[n3d] = get_next_token_int(&token) - 1;
            nmat3d[sm->elem3d_physics_mat[n3d]]++;
            n3d++;
        } else if ( (strcmp(token,"TRI")  == 0) || (strcmp(token,"QUAD") == 0) ) {
            sm->elem2d_physics_mat[n2d] = get_next_token_int(&token) - 1;
            nmat2d[sm->elem2d_physics_mat[n2d]]++;
            n2d++;
        } else if ( (strcmp(token,"SEG") == 0) ) {
            sm->elem1d_physics_mat[n1d] = get_next_token_int(&token) - 1;  
            nmat1d[sm->elem1d_physics_mat[n1d]]++;
            n1d++;
        }
    }
    assert(n1d == sm->grid->nelems1d); 
    assert(n2d == sm->grid->nelems2d);
    assert(n3d == sm->grid->nelems3d);
    fclose(fp);
    if (DEBUG) {
        printf("Physics Coverage File Info:\n");
        printf("-- n1d: %d || nelems1d: %d\n",n1d,sm->grid->nelems1d);
        printf("-- n2d: %d || nelems2d: %d\n",n2d,sm->grid->nelems2d);
        printf("-- n3d: %d || nelems3d: %d\n",n3d,sm->grid->nelems3d);
        for (i=0; i<30; i++) {if (nmat3d[i] > 0) printf("-- 3D physics material[%d] used %d times\n",i,nmat3d[i]);}
        for (i=0; i<30; i++) {if (nmat2d[i] > 0) printf("-- 2D physics material[%d] used %d times\n",i,nmat2d[i]);}
        for (i=0; i<30; i++) {if (nmat1d[i] > 0) printf("-- 1D physics material[%d] used %d times\n",i,nmat1d[i]);}
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Create Node Materials Based by Pointing to Highest DOF SMAT_PHYSICS
    //++++++++++++++++++++++++++++++++++++++++++++++

    if (DEBUG) {
        printf(">creating mat_physics_node\n");
    }
    sm->mat_physics_node = (SMAT_PHYSICS **) tl_alloc(sizeof(SMAT_PHYSICS *),grid->nnodes);
    int node_nvars[grid->nnodes]; sarray_init_int(node_nvars,grid->nnodes);
    
    int ie,imat,nd;
    for (ie=0; ie<grid->nelems1d; ie++) {
        if (grid->elem1d[ie].bflag != BODY) continue;
        imat = sm->elem1d_physics_mat[ie];
        for (i=0; i<grid->elem1d[ie].nnodes; i++) {
            nd = grid->elem1d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem1d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems2d; ie++) {
        if (grid->elem2d[ie].bflag != BODY) continue;
        imat = sm->elem2d_physics_mat[ie];
        for (i=0; i<grid->elem2d[ie].nnodes; i++) {
            nd = grid->elem2d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem2d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    for (ie=0; ie<grid->nelems3d; ie++) {
        imat = sm->elem3d_physics_mat[ie];
        for (i=0; i<grid->elem3d[ie].nnodes; i++) {
            nd = grid->elem3d[ie].nodes[i];
            if (node_nvars[nd] < sm->mat_physics_elem[imat].ivar_pos.n) {
                node_nvars[nd] = sm->mat_physics_elem[imat].ivar_pos.n;
                sm->mat_physics_node[nd] =  &(sm->mat_physics_elem[sm->elem3d_physics_mat[ie]]); // point to the elemental material
            }
        }
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Find all independent variables used in the SuperModel
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">creating ivar_pos\n");
    }
    int FLAGS[MAX_TRNS_VARS + MAX_VARS]; sarray_init_int(FLAGS,MAX_TRNS_VARS + MAX_VARS);
    sivar_position_init(&(sm->ivar_pos));
    smat_physics_position_flag(sm->mat_physics_node,grid->nnodes,FLAGS); 
    sivar_position_map(&(sm->ivar_pos),FLAGS);
    //printf("sm->ivar_pos.n: %d\n",sm->ivar_pos.n);
    //printf("sm->ivar_pos.ntrns: %d\n",sm->ivar_pos.ntrns);
    if (DEBUG) {
        printf("Variable Position Info:\n");
        sivar_position_printScreen(&sm->ivar_pos);
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Use this to update ivar pos in mat physics
    //++++++++++++++++++++++++++++++++++++++++++++++
    smat_physics_update_array(sm->mat_physics_elem,sm->nmat_physics,&(sm->ivar_pos));

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate the map from the independent variables 
    // to their position in the solution vector
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf(">allocating solution variables\n");
    }
    // problem here
    SIVAR_POSITION ip[grid->nnodes];
    for (i=0; i<grid->nnodes; i++) {ip[i] = sm->mat_physics_node[i]->ivar_pos;}
    sm->ndofs = sivar_position_build_dof_map(&sm->ivar_pos,grid->nnodes,ip,&sm->ivars);
    if (DEBUG) {
        printf(">total ndofs: %d\n",sm->ndofs);
    }
    sm->sol       = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_old   = (double *) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol_older = (double *) tl_alloc(sizeof(double), sm->ndofs);

    sarray_init_dbl(sm->sol,sm->ndofs);
    sarray_init_dbl(sm->sol_old,sm->ndofs);
    sarray_init_dbl(sm->sol_older,sm->ndofs);
    
}

