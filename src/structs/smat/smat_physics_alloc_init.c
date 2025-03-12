/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smat_physics_alloc_init The file is for allocating an initializing an SMAT_PHYSICS struct */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
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

void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, char **codes) {
    assert(nmat >= 0);
    if (nmat == 0){
        mat_physics = NULL;
        return;
    } else {
        (*mat_physics) = (SMAT_PHYSICS *) tl_alloc(sizeof(SMAT_PHYSICS), nmat);
        SMAT_PHYSICS *mat;  // alias
        for (int imat=0; imat<nmat; imat++) {
            mat = &(*mat_physics)[imat];
            mat->id = imat;
            smat_physics_alloc_init(mat, codes[imat]);
        }
        return;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 * @param[in]    codes       (char [][10]) an array of strings that encode the materials physics
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_alloc_init(SMAT_PHYSICS *mat, char *codes) {

    strcpy(mat->code,codes);
    if (strcmp(mat->code,"0000") == 0) {
        printf("ERROR: No physics given on material %d\n",mat->id);
        tl_error("");
    }
    sivar_position_init(&mat->ivar_pos);

    int itrns;

    // get number of transport contistuents on this element
    mat->ivar_pos.ntrns = mat->code[3] - '0'; // only 4th column for now [0-9] transports
    mat->TRANSPORT = NULL;
    //printf("mat->ivar_pos.ntrns: %d\n",mat->ivar_pos.ntrns);
    if (mat->ivar_pos.ntrns > 0) {
        mat->TRANSPORT = (bool *) tl_alloc(sizeof(bool), mat->ivar_pos.ntrns);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    // initialize all models to off
    //++++++++++++++++++++++++++++++++++++++++++++++
    // independent variables
    mat->SW_FLOW    = false;
    mat->SW1_FLOW   = false;
    mat->SW2_FLOW   = false;
    mat->SW3_FLOW   = false;
    mat->NS_FLOW    = false;
    mat->NS3_FLOW   = false;
    mat->NS3_SPLIT  = false;
    mat->GW_FLOW    = false;
    mat->DW_FLOW    = false;
    mat->WVEL_SPLIT = false;
    mat->PRESSURE   = false;
    mat->POISSON    = false;
    for (itrns=0; itrns<mat->ivar_pos.ntrns; itrns++) {
        mat->TRANSPORT[itrns] = false;
    }
        // dependent variables
    mat->VORTICITY = false;
    mat->SEDIMENT  = false;
    mat->SEDLIB    = false;
    mat->ICM       = false;
    mat->NSM       = false;
    mat->WAVE      = false;
    mat->WIND      = false;

        //++++++++++++++++++++++++++++++++++++++++++++++
        // load the variables into the array
        //++++++++++++++++++++++++++++++++++++++++++++++
    mat->ivar_pos.n = 0;
    mat->nSubmodels = 0;

    // surface water
    if (mat->code[0] == '1') {
        mat->SW_FLOW = true;
        mat->SW1_FLOW = true;
        mat->ivar_pos.var[_H]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 1D SW
        mat->ivar_pos.var[_UDA] = mat->ivar_pos.n; mat->ivar_pos.n++; // 1D SW
        mat->nSubmodels++;
    } else if (mat->code[0] == '2') {
        mat->SW_FLOW = true;
        mat->SW2_FLOW = true;
        mat->ivar_pos.var[_H]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 2D SW
        mat->ivar_pos.var[_UDA] = mat->ivar_pos.n; mat->ivar_pos.n++; // 2D SW
        mat->ivar_pos.var[_VDA] = mat->ivar_pos.n; mat->ivar_pos.n++; // 2D SW
        mat->nSubmodels++;
    } else if (mat->code[0] == '3') {
        mat->SW_FLOW = true;
        mat->SW3_FLOW = true;
        mat->ivar_pos.var[_DPL] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D SW
        mat->ivar_pos.var[_U]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D SW
        mat->ivar_pos.var[_V]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D SW
        mat->nSubmodels++;
    } else if (mat->code[0] == '4') {
        mat->NS_FLOW = true;
        mat->NS3_FLOW = true;
        mat->ivar_pos.var[_U]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS
        mat->ivar_pos.var[_V]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS
        mat->ivar_pos.var[_W]   = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS
        mat->ivar_pos.var[_PRS] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS
        mat->nSubmodels++;
    } else if (mat->code[0] == '5') {
        mat->NS_FLOW = true;
        mat->NS3_SPLIT = true;
        mat->ivar_pos.var[_U] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS SPLIT
        mat->ivar_pos.var[_V] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS SPLIT
        mat->ivar_pos.var[_W] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS SPLIT
        mat->nSubmodels++;
    } else if (mat->code[0] == '6') {
        mat->DW_FLOW = true;
        mat->ivar_pos.var[_H] = mat->ivar_pos.n; mat->ivar_pos.n++;   // 2D OF
        mat->nSubmodels++;
    } else if (mat->code[0] == '7') {
        mat->WVEL_SPLIT = true;
        mat->ivar_pos.var[_W] = mat->ivar_pos.n; mat->ivar_pos.n++;   // 3D SW - SPLIT
        mat->nSubmodels++;
    } else if (mat->code[0] == '8') {
        mat->PRESSURE = true;
        mat->ivar_pos.var[_PRS] = mat->ivar_pos.n; mat->ivar_pos.n++; // 3D NS - SPLIT
        mat->nSubmodels++;
    }else if (mat->code[0] == '9'){
        //just for testing
        mat->POISSON = true;
        mat->ivar_pos.var[_H] = mat->ivar_pos.n; mat->ivar_pos.n++; // 2D Poisson
        mat->nSubmodels++;
    }
    
    // ground water
    if (mat->code[1] == '1') {
        mat->GW_FLOW = true;
        // only add depth if surface water has not
        if (mat->ivar_pos.var[_H]  == UNSET_INT) {
            mat->ivar_pos.var[_H] = mat->ivar_pos.n; mat->ivar_pos.n++; // 1D/2D/3D GROUNDWATER
        }
        mat->nSubmodels++;
    }
    
    // transport
    mat->ntrns = 0;
    for (itrns=0; itrns<mat->ivar_pos.ntrns; itrns++) {
        mat->TRANSPORT[itrns] = true;
        mat->ivar_pos.var[N_IVARS + itrns] = mat->ivar_pos.n;
        mat->ntrns++;
        mat->ivar_pos.n++; // 1D/2D/3D TRANSPORT
        mat->nSubmodels++;
    }

    assert(mat->nSubmodels > 0);
    
    //Mark adding n, this is same as before
    //can consolidate later
    mat->n = mat->ivar_pos.n;
    mat->ivar_loc = (int *) tl_alloc(sizeof(int), mat->n);
    sarray_init_value_int(mat->ivar_loc, mat->n, UNSET_INT);
    int nSubMod_nvar[mat->nSubmodels];
    sarray_init_int(nSubMod_nvar,mat->nSubmodels);



    int k=0;
    if (mat->SW1_FLOW)   {nSubMod_nvar[k]=2; k++;} // 1
    if (mat->SW2_FLOW)   {nSubMod_nvar[k]=3; k++;} // 2
    if (mat->SW3_FLOW)   {nSubMod_nvar[k]=3; k++;} // 3
    if (mat->NS3_FLOW)   {nSubMod_nvar[k]=4; k++;} // 4
    if (mat->NS3_SPLIT)  {nSubMod_nvar[k]=3; k++;} // 5
    if (mat->DW_FLOW)    {nSubMod_nvar[k]=1; k++;} // 6
    if (mat->WVEL_SPLIT) {nSubMod_nvar[k]=1; k++;} // 7
    if (mat->PRESSURE)   {nSubMod_nvar[k]=1; k++;} // 8
    if (mat->POISSON)    {nSubMod_nvar[k]=1; k++; printf("Set POISONN \n\n\n");} // 9
    if (mat->GW_FLOW)    {nSubMod_nvar[k]=1; k++;}
    for (itrns=0; itrns<mat->ivar_pos.ntrns; itrns++) {
        if (mat->TRANSPORT[itrns]) {nSubMod_nvar[k]=1; k++;}    
    }

    int tmp[mat->nSubmodels], tmp2[mat->nSubmodels]; // Just to allocate
    sarray_init_value_int(tmp,mat->nSubmodels,UNSET_INT);
    sarray_init_value_int(tmp2,mat->nSubmodels,UNSET_INT);
    smodel_alloc_init_array(&(mat->model), tmp, tmp2, mat->nSubmodels, nSubMod_nvar); 

    //Remainder of models will be filled in in smat_physics_update_array.c
}

