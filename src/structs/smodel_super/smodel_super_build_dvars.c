/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super_build_dvars.c This file collects methods of the SCOVERAGE structure for coverages   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void add_dvar(int just_count,SDVAR_POSITION *pos,int ID,int *n);
void load_dvars(SMODEL_SUPER *sm,int *nvar_node, int *nvar_elem_dbl, int *nvar_elem_int, bool just_count);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Fills out the DVARS structure of a SuperModel instance for storing all require dependent variables
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SMODEL_SUPER *)  a pointer to a superModel

 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_build_dvars(SMODEL_SUPER *sm) {
    
    int nnodes = sm->grid->nnodes;
    int nnode_dvar = sm->grid->nnodes;   // FOR NOW!
    int nelem_dvar = sm->grid->nelems2d; // FOR NOW!
    
    int nvar_node = 0, nvar_elem_dbl = 0, nvar_elem_int = 0;
    load_dvars(sm,&nvar_node,&nvar_elem_dbl,&nvar_elem_int,true);
    sdvar_alloc_init(&sm->dvars,nnodes,nvar_node,nvar_elem_dbl,nvar_elem_int,nnode_dvar,nelem_dvar);
    load_dvars(sm,&nvar_node,&nvar_elem_dbl,&nvar_elem_int,false);
    
    if (DEBUG) {
        printf("nvar_node: %d || nvar_elem_dbl: %d || nvar_elem_int: %d\n",nvar_node,nvar_elem_dbl,nvar_elem_int);
        if (nvar_node > 0)     {sdvar_position_printScreen(&sm->dvars.sdvar_pos_node);}
        if (nvar_elem_int > 0) {sdvar_position_printScreen(&sm->dvars.sdvar_pos_elem_int);}
        if (nvar_elem_dbl > 0) {sdvar_position_printScreen(&sm->dvars.sdvar_pos_elem_dbl);}
    }
    
}

void add_dvar(int just_count,SDVAR_POSITION *pos,int ID,int *n) {
    if (!just_count) {pos->var[ID] = *n;}
    (*n)++;
    pos->n = *n;
}

void load_dvars(SMODEL_SUPER *sm,int *nvar_node, int *nvar_elem_dbl, int *nvar_elem_int, bool just_count) {
    
    // aliases
    SDVAR_POSITION *npos  = &(sm->dvars.sdvar_pos_node);
    SDVAR_POSITION *eposi = &(sm->dvars.sdvar_pos_elem_int);
    SDVAR_POSITION *eposd = &(sm->dvars.sdvar_pos_elem_dbl);
    
    *nvar_node = 0; *nvar_elem_dbl = 0; *nvar_elem_int = 0;
    
    // use flags to decide which dependent variables are needed
    if (sm->flags.model[adh_def._SW1] == 1 ||
        sm->flags.model[adh_def._SW2] == 1 ||
        sm->flags.model[adh_def._SW3] == 1) {
        add_dvar(just_count,npos,adh_def._DPRS,nvar_node);
        add_dvar(just_count,npos,adh_def._DENSITY,nvar_node);
        add_dvar(just_count,npos,adh_def._ERROR,nvar_node);
        
        add_dvar(just_count,eposi,adh_def._WDFLAG,nvar_elem_int);
        add_dvar(just_count,eposd,adh_def._ERROR,nvar_elem_dbl);
    }
    if (sm->flags.WIND_STRESS) {
        add_dvar(just_count,npos,adh_def._WIND_SX,nvar_node);
        add_dvar(just_count,npos,adh_def._WIND_SY,nvar_node);
    }
    if (sm->flags.WAVE_STRESS) {
        add_dvar(just_count,npos,adh_def._WAVE_SX,nvar_node);
        add_dvar(just_count,npos,adh_def._WAVE_SY,nvar_node);
    }
    if (sm->flags.WAVE_RADS) {
        add_dvar(just_count,npos,adh_def._WAVE_XX,nvar_node);
        add_dvar(just_count,npos,adh_def._WAVE_XY,nvar_node);
        add_dvar(just_count,npos,adh_def._WAVE_YX,nvar_node);
        add_dvar(just_count,npos,adh_def._WAVE_YY,nvar_node);
    }
     
    
}

