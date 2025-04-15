/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sadh\_def.c Collects routine for AdH models and variables  */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Initializes a SADH\_DEF structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] def (SADH_DEF *)  a pointer to a SADH_DEF structure that stores AdH model and variable dictionaries
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sadh_def_init(void) {
    int i=0, k=0, imod=0;
    
    // ++++++++++++++++++++++++++++++++++
    // ADH INDEPENDENT VARIABLES
    // ++++++++++++++++++++++++++++++++++
    k=0;
    adh_def.n_ivars = 10 + MAX_TRNS;
    adh_def.ivar = (SVAR *) tl_alloc(sizeof(SVAR),adh_def.n_ivars);
    adh_def._H     = k; strcpy(adh_def.ivar[k].name,"DEPTH");       strcpy(adh_def.ivar[k].subname,"");  adh_def.ivar[k].dim = 1; k++;
    adh_def._U     = k; strcpy(adh_def.ivar[k].name,"VELOCITY");    strcpy(adh_def.ivar[k].subname,"X"); adh_def.ivar[k].dim = 3; k++;
    adh_def._V     = k; strcpy(adh_def.ivar[k].name,"VELOCITY");    strcpy(adh_def.ivar[k].subname,"Y"); adh_def.ivar[k].dim = 3; k++;
    adh_def._W     = k; strcpy(adh_def.ivar[k].name,"VELOCITY");    strcpy(adh_def.ivar[k].subname,"Z"); adh_def.ivar[k].dim = 3; k++;
    adh_def._UDA   = k; strcpy(adh_def.ivar[k].name,"DA_VELOCITY"); strcpy(adh_def.ivar[k].subname,"X"); adh_def.ivar[k].dim = 2; k++;
    adh_def._VDA   = k; strcpy(adh_def.ivar[k].name,"DA_VELOCITY"); strcpy(adh_def.ivar[k].subname,"Y"); adh_def.ivar[k].dim = 2; k++;
    adh_def._DPL   = k; strcpy(adh_def.ivar[k].name,"DISPLACEMENT");strcpy(adh_def.ivar[k].subname,"");  adh_def.ivar[k].dim = 1; k++;
    adh_def._PRS   = k; strcpy(adh_def.ivar[k].name,"PRESSURE");    strcpy(adh_def.ivar[k].subname,"");  adh_def.ivar[k].dim = 1; k++;
    adh_def._HEAT  = k; strcpy(adh_def.ivar[k].name,"HEAT");        strcpy(adh_def.ivar[k].subname,"");  adh_def.ivar[k].dim = 1; k++;
    adh_def._SAL   = k; strcpy(adh_def.ivar[k].name,"SALINITY");    strcpy(adh_def.ivar[k].subname,"");  adh_def.ivar[k].dim = 1; k++;
    for (i=0; i<MAX_TRNS; i++) {
        adh_def._CON[i]  = k; sprintf(adh_def.ivar[k].name, "CON%d", i); strcpy(adh_def.ivar[k].subname,""); adh_def.ivar[k].dim = 1; k++;
    }
    assert(k == adh_def.n_ivars);
    
    // ++++++++++++++++++++++++++++++++++
    // ADH INDEPENDENT VARIABLES
    // ++++++++++++++++++++++++++++++++++
    k=0;
    adh_def.n_dvars = 19;
    adh_def.dvar = (SVAR *) tl_alloc(sizeof(SVAR),adh_def.n_dvars);
    adh_def._WIND_SX    = k; strcpy(adh_def.dvar[k].name,"WIND_STRESS");       strcpy(adh_def.dvar[k].subname,"X");   adh_def.dvar[k].dim = 2; k++;
    adh_def._WIND_SY    = k; strcpy(adh_def.dvar[k].name,"WIND_STRESS");       strcpy(adh_def.dvar[k].subname,"Y");   adh_def.dvar[k].dim = 2; k++;
    adh_def._WAVE_SX    = k; strcpy(adh_def.dvar[k].name,"WAVE_STRESS");       strcpy(adh_def.dvar[k].subname,"X");   adh_def.dvar[k].dim = 2; k++;
    adh_def._WAVE_SY    = k; strcpy(adh_def.dvar[k].name,"WAVE_STRESS");       strcpy(adh_def.dvar[k].subname,"Y");   adh_def.dvar[k].dim = 2; k++;
    adh_def._WAVE_XX    = k; strcpy(adh_def.dvar[k].name,"WAVE_RADS");         strcpy(adh_def.dvar[k].subname,"XX");  adh_def.dvar[k].dim = 4; k++;
    adh_def._WAVE_XY    = k; strcpy(adh_def.dvar[k].name,"WAVE_RADS");         strcpy(adh_def.dvar[k].subname,"XY");  adh_def.dvar[k].dim = 4; k++;
    adh_def._WAVE_YX    = k; strcpy(adh_def.dvar[k].name,"WAVE_RADS");         strcpy(adh_def.dvar[k].subname,"YX");  adh_def.dvar[k].dim = 4; k++;
    adh_def._WAVE_YY    = k; strcpy(adh_def.dvar[k].name,"WAVE_RADS");         strcpy(adh_def.dvar[k].subname,"YY");  adh_def.dvar[k].dim = 4; k++;
    adh_def._DENSITY    = k; strcpy(adh_def.dvar[k].name,"DENSITY");           strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._BED_DPL    = k; strcpy(adh_def.dvar[k].name,"BED_DISPLACEMENT");  strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._BED_ELE    = k; strcpy(adh_def.dvar[k].name,"BED_ELEVATION");     strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._GS         = k; strcpy(adh_def.dvar[k].name,"GRID_SPEED");        strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._GS_OLD     = k; strcpy(adh_def.dvar[k].name,"GRID_SPEED_OLD");    strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._DPRS       = k; strcpy(adh_def.dvar[k].name,"PRESSURE");          strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._DPRS_PLUS  = k; strcpy(adh_def.dvar[k].name,"PRS_PLUS");          strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._DPRS_MINUS = k; strcpy(adh_def.dvar[k].name,"PRS_MINUS");         strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._VNODE_FLUX = k; strcpy(adh_def.dvar[k].name,"VERTICLE_NODE_FLUX");strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._WDFLAG     = k; strcpy(adh_def.dvar[k].name,"WET_DRY_FLAG");      strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    adh_def._ERROR      = k; strcpy(adh_def.dvar[k].name,"ERROR");             strcpy(adh_def.dvar[k].subname,"");    adh_def.dvar[k].dim = 1; k++;
    assert(k == adh_def.n_dvars);
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ADH MODEL DEFINITIONS
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // CJT - note each transport is a separate model
    adh_def.nmodels = 11 + MAX_TRNS * 3;
    adh_def.model = (SADH_MODEL_DEF *) tl_alloc(sizeof(SADH_MODEL_DEF), adh_def.nmodels);
    for (imod=0; imod<adh_def.nmodels; imod++) {
        strcpy(adh_def.model[imod].name,"");
        strcpy(adh_def.model[imod].subname,"");
        sarray_init_int(adh_def.model[imod].var,MAX_IVAR_DIM);
    }
    
    imod=0;
    
    //+++++++++++++++++++++++++++++++++++
    // 1D SW ++++++++++++++++++++++++++++
    adh_def._SW1 = imod;
    strcpy(adh_def.model[imod].name,"SW");
    strcpy(adh_def.model[imod].subname,"SW1");
    adh_def.model[imod].nivars = 2;
    adh_def.model[imod].var[0]   = adh_def._H;
    adh_def.model[imod].var[1]   = adh_def._UDA;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D SW ++++++++++++++++++++++++++++
    adh_def._SW2 = imod;
    strcpy(adh_def.model[imod].name,"SW");
    strcpy(adh_def.model[imod].subname,"SW2");
    adh_def.model[imod].nivars = 3;
    adh_def.model[imod].var[0]   = adh_def._H;
    adh_def.model[imod].var[1]   = adh_def._UDA;
    adh_def.model[imod].var[2]   = adh_def._VDA;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 3D SW ++++++++++++++++++++++++++++
    adh_def._SW3 = imod;
    strcpy(adh_def.model[imod].name,"SW");
    strcpy(adh_def.model[imod].subname,"SW3");
    adh_def.model[imod].nivars = 3;
    adh_def.model[imod].var[0]   = adh_def._H;
    adh_def.model[imod].var[1]   = adh_def._U;
    adh_def.model[imod].var[2]   = adh_def._V;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 3D NS ++++++++++++++++++++++++++++
    adh_def._NS3 = imod;
    strcpy(adh_def.model[imod].name,"NS");
    strcpy(adh_def.model[imod].subname,"NS3");
    adh_def.model[imod].nivars = 4;
    adh_def.model[imod].var[0]   = adh_def._U;
    adh_def.model[imod].var[1]   = adh_def._V;
    adh_def.model[imod].var[2]   = adh_def._W;
    adh_def.model[imod].var[3]   = adh_def._PRS;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 3D NS SPLIT ++++++++++++++++++++++
    adh_def._NS3SPLIT = imod;
    strcpy(adh_def.model[imod].name,"NS");
    strcpy(adh_def.model[imod].subname,"NS3SPLIT");
    adh_def.model[imod].nivars = 4;
    adh_def.model[imod].var[0]   = adh_def._U;
    adh_def.model[imod].var[1]   = adh_def._V;
    adh_def.model[imod].var[2]   = adh_def._W;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D DW ++++++++++++++++++++++++++++
    adh_def._DW2D = imod;
    strcpy(adh_def.model[imod].name,"DW");
    strcpy(adh_def.model[imod].subname,"DW2");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._H;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D WVEL ++++++++++++++++++++++++++
    adh_def._WVEL = imod;
    strcpy(adh_def.model[imod].name,"WVEL");
    strcpy(adh_def.model[imod].subname,"WVEL");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._W;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D PRS ++++++++++++++++++++++++++
    adh_def._PRS = imod;
    strcpy(adh_def.model[imod].name,"PRESSURE");
    strcpy(adh_def.model[imod].subname,"PRESSURE");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._PRS;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D POISSON ++++++++++++++++++++++
    adh_def._POISSON = imod;
    strcpy(adh_def.model[imod].name,"POISSON");
    strcpy(adh_def.model[imod].subname,"POISSON2D");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._H;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D HEAT +++++++++++++++++++++++++
    adh_def._HEAT2D = imod;
    strcpy(adh_def.model[imod].name,"HEAT");
    strcpy(adh_def.model[imod].subname,"HEAT2D");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._H;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // 2D GW3D +++++++++++++++++++++++++
    adh_def._GW3D = imod;
    strcpy(adh_def.model[imod].name,"GW");
    strcpy(adh_def.model[imod].subname,"GW3D");
    adh_def.model[imod].nivars = 1;
    adh_def.model[imod].var[0]   = adh_def._H;
    imod++;
    
    //+++++++++++++++++++++++++++++++++++
    // TRANSPORT +++++++++++++++++++++++
    for (i=0; i<MAX_TRNS; i++) {
        adh_def._TRNS1D[i] = imod;
        strcpy(adh_def.model[imod].name,"TRANSPORT");
        strcpy(adh_def.model[imod].subname,"TRNS1D");
        adh_def.model[imod].nivars = 1;
        adh_def.model[imod].var[0]   = adh_def._CON[i];
        imod++;
        
        adh_def._TRNS2D[i] = imod;
        strcpy(adh_def.model[imod].name,"TRANSPORT");
        strcpy(adh_def.model[imod].subname,"TRNS2D");
        adh_def.model[imod].nivars = 1;
        adh_def.model[imod].var[0]   = adh_def._CON[i];
        imod++;
        
        adh_def._TRNS3D[i] = imod;
        strcpy(adh_def.model[imod].name,"TRANSPORT");
        strcpy(adh_def.model[imod].subname,"TRNS3D");
        adh_def.model[imod].nivars = 1;
        adh_def.model[imod].var[0]   = adh_def._CON[i];
        imod++;
    }
    printf("imod: %d || adh_def.nmodels: %d \n",imod,adh_def.nmodels);
    assert(imod == adh_def.nmodels);
}

void sadh_def_sivar_printScreen(SVAR v) {
    printf(">>>>> IVAR || name: %s - %s || dim: %d\n",v.name,v.subname,v.dim);
}

void sadh_def_model_printScreen(SADH_MODEL_DEF *m) {
    printf(">>> SMODEL %d || name: %s - %s || nivars: %d\n",m->id,m->name,m->subname,m->nivars);
    for (int i=0; i<m->nivars; i++) {
        sadh_def_sivar_printScreen(adh_def.ivar[m->var[i]]);
    }
}

residual_ptr select_resid_func(char *model, int *imod, int trns_id) {
    if (strcmp(model,"SW2") == 0) {
        *imod = adh_def._SW2;
        return fe_sw2_body_resid;
    } else if (strcmp(model,"DW") == 0) {
        *imod = adh_def._DW2D;
        return NULL;
    } else if (strcmp(model,"POISSON2D") == 0) {
        *imod = adh_def._POISSON;
        return poisson_residual;
    } else if (strcmp(model,"TRNS") == 0) {
        *imod = adh_def._TRNS2D[trns_id];
        return NULL;
    } else if (strcmp(model, "HEAT2D") == 0){
        *imod = adh_def._HEAT2D;
        return heat_residual;
    } else {
        *imod = UNSET_INT;
        return NULL;
    }
}

residual_ptr select_bresid_func(char *model, char *bc_type, char *bc_var, int *imod, int trns_id) {
        *imod = UNSET_INT;
        if (strcmp(model,"SW2") == 0) {
            if (strcmp(bc_type,"NB") == 0) {
                if (strcmp(bc_var,"VEL") == 0) {
                    *imod = adh_def._SW2;
                    return fe_sw2_bc_vel;
                }
            } else if (strcmp(bc_type,"DEPTH") == 0) {
                *imod = adh_def._SW2;
                return fe_sw2_bc_h;
            } else {
                return NULL;
            }
    
        } else if (strcmp(model,"POISSON2D") == 0) {
            *imod = adh_def._POISSON;
            return NULL;
    
        } else {
            return NULL;
        }
}

init_func_ptr select_init_func(char *model, int trns_id) {
    if (strcmp(model,"SW2") == 0) {
        return fe_sw2_init;
    } else if (strcmp(model,"DW") == 0) {
        return NULL;
    } else if (strcmp(model,"TRNS") == 0) {
        return NULL;
    } else if (strcmp(model,"POISSON2D") == 0) {
        return NULL;
    } else {
        return NULL;
    }

}
