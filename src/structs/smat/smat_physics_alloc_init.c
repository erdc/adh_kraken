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

void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, char codes[][10]) {
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
void smat_physics_alloc_init(SMAT_PHYSICS *mat, char codes[10]) {

    strcpy(mat->code,codes);
    if (strcmp(mat->code,"0000") == 0) {
        printf("ERROR: No physics given on material %d\n",mat->id);
        tl_error("");
    }
    sivar_position_init(&mat->ivars);

    int itrns;

    // get number of transport contistuents on this element
    mat->ivars.ntrns = mat->code[3] - '0'; // only 4th column for now [0-9] transports
    //printf("mat->ivars.ntrns: %d\n",mat->ivars.ntrns);
    if (mat->ivars.ntrns > 0) {
        mat->TRANSPORT = (bool *) tl_alloc(sizeof(bool), mat->ivars.ntrns);
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
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
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
    mat->ivars.n = 0;
    mat->nSubmodels = 0;

        // surface water
    if (mat->code[0] == '1') {
        mat->SW_FLOW = true;
        mat->SW1_FLOW = true;
            mat->ivars.h   = mat->ivars.n; mat->ivars.n++; // 1D SW
            mat->ivars.uda = mat->ivars.n; mat->ivars.n++; // 1D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '2') {
            mat->SW_FLOW = true;
            mat->SW2_FLOW = true;
            mat->ivars.h   = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->ivars.uda = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->ivars.vda = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '3') {
            mat->SW_FLOW = true;
            mat->SW3_FLOW = true;
            mat->ivars.dpl = mat->ivars.n; mat->ivars.n++; // 3D SW 
            mat->ivars.u   = mat->ivars.n; mat->ivars.n++; // 3D SW
            mat->ivars.v   = mat->ivars.n; mat->ivars.n++; // 3D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '4') {
            mat->NS_FLOW = true;
            mat->NS3_FLOW = true;
            mat->ivars.u   = mat->ivars.n;   mat->ivars.n++; // 3D NS 
            mat->ivars.v   = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->ivars.w   = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->ivars.prs = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->nSubmodels++;
        } else if (mat->code[0] == '5') {
            mat->NS_FLOW = true;
            mat->NS3_SPLIT = true;
            mat->ivars.u = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->ivars.v = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->ivars.w = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->nSubmodels++;
        } else if (mat->code[0] == '6') {
            mat->DW_FLOW = true;
            mat->ivars.h = mat->ivars.n; mat->ivars.n++;   // 2D OF
            mat->nSubmodels++;
        } else if (mat->code[0] == '7') {
            mat->WVEL_SPLIT = true;
            mat->ivars.w = mat->ivars.n; mat->ivars.n++;   // 3D SW - SPLIT
            mat->nSubmodels++;
        } else if (mat->code[0] == '8') {
            mat->PRESSURE = true;
            mat->ivars.prs = mat->ivars.n; mat->ivars.n++; // 3D NS - SPLIT
            mat->nSubmodels++;
        }

        // ground water
        if (mat->code[1] == '1') {
            mat->GW_FLOW = true;
            // only add depth if surface water has not
            if (mat->ivars.h == UNSET_INT) {
                mat->ivars.h = mat->ivars.n; mat->ivars.n++; // 1D/2D/3D GROUNDWATER
            }
            mat->nSubmodels++;
        }

        // transport
        mat->ntrns = 0;
        for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
            mat->TRANSPORT[itrns] = true;
            mat->ntrns++;
            mat->ivars.con[itrns] = mat->ivars.n; mat->ivars.n++; // 1D/2D/3D TRANSPORT
            mat->nSubmodels++;
        }
    assert(mat->nSubmodels > 0);
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
    if (mat->GW_FLOW)    {nSubMod_nvar[k]=1; k++;}
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
        if (mat->TRANSPORT[itrns]) {nSubMod_nvar[k]=1; k++;}    
    }

    int tmp[mat->nSubmodels], tmp2[mat->nSubmodels]; // Just to allocate
    sarray_init_value_int(tmp,mat->nSubmodels,UNSET_INT);
    sarray_init_value_int(tmp2,mat->nSubmodels,UNSET_INT);
    smodel_alloc_init_array(&(mat->model), tmp, tmp2, mat->nSubmodels, nSubMod_nvar); 

    // fill in model info
    int isubModel = 0;
    if (mat->SW1_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.uda;
        mat->model[isubModel].physics = SW1D_; // for body residuals
        mat->model[isubModel].physics_init = SW1D_; 
        isubModel++;
    }
    if (mat->SW2_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.uda;
        mat->model[isubModel].physics_vars[2] = mat->ivars.vda;
        mat->model[isubModel].physics = SW2D_;
        mat->model[isubModel].physics_init = SW2D_;
        isubModel++;
    }
    if (mat->SW3_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.u;
        mat->model[isubModel].physics_vars[2] = mat->ivars.v;
        mat->model[isubModel].physics = SW3D_;
        mat->model[isubModel].physics_init = SW3D_;
        isubModel++;
    }
    if (mat->NS3_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.u;
        mat->model[isubModel].physics_vars[1] = mat->ivars.v;
        mat->model[isubModel].physics_vars[2] = mat->ivars.w;
        mat->model[isubModel].physics_vars[3] = mat->ivars.prs;
        mat->model[isubModel].physics = NS3D_;
        mat->model[isubModel].physics_init = NS3D_;
        isubModel++;
    }
    if (mat->NS3_SPLIT) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.u;
        mat->model[isubModel].physics_vars[1] = mat->ivars.v;
        mat->model[isubModel].physics_vars[2] = mat->ivars.w;
        mat->model[isubModel].physics = NS3DSPLIT_;
        mat->model[isubModel].physics_init = NS3DSPLIT_;
        isubModel++;
    }
    if (mat->DW_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics = DW2D_;
        mat->model[isubModel].physics_init = DW2D_;
        isubModel++;
    }
    if (mat->WVEL_SPLIT) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.w;
        mat->model[isubModel].physics = WVEL_;
        mat->model[isubModel].physics_init = WVEL_;
        isubModel++;
    }
    if (mat->PRESSURE) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.prs;
        mat->model[isubModel].physics = PRS_;
        mat->model[isubModel].physics_init = PRS_;
        isubModel++;
    }
    if (mat->GW_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics = GW3D_;
        mat->model[isubModel].physics_init = GW3D_;
        isubModel++;
    }
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
        if (mat->TRANSPORT[itrns]) {
            mat->model[isubModel].physics_vars[0] = mat->ivars.con[itrns];
            mat->model[isubModel].physics = TRNS_;
            mat->model[isubModel].physics_init = TRNS_;
            isubModel++;
        }  
    }
    assert(isubModel == mat->nSubmodels);
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
 * @param[in]    nmat        (int) number of materials to allocate
 * @param[in]    codes       (char [][10]) an array of strings that encode the materials physics
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_physics_alloc_init_ptr_array(SMAT_PHYSICS **mat_physics, int nmat, char **codes) {
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
            smat_physics_alloc_init_ptr(mat, codes[imat]);
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
void smat_physics_alloc_init_ptr(SMAT_PHYSICS *mat, char *codes) {

    strcpy(mat->code,codes);
    if (strcmp(mat->code,"0000") == 0) {
        printf("ERROR: No physics given on material %d\n",mat->id);
        tl_error("");
    }
    sivar_position_init(&mat->ivars);

    int itrns;

    // get number of transport contistuents on this element
    mat->ivars.ntrns = mat->code[3] - '0'; // only 4th column for now [0-9] transports
    //printf("mat->ivars.ntrns: %d\n",mat->ivars.ntrns);
    if (mat->ivars.ntrns > 0) {
        mat->TRANSPORT = (bool *) tl_alloc(sizeof(bool), mat->ivars.ntrns);
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
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
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
    mat->ivars.n = 0;
    mat->nSubmodels = 0;

        // surface water
    if (mat->code[0] == '1') {
        mat->SW_FLOW = true;
        mat->SW1_FLOW = true;
            mat->ivars.h   = mat->ivars.n; mat->ivars.n++; // 1D SW
            mat->ivars.uda = mat->ivars.n; mat->ivars.n++; // 1D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '2') {
            mat->SW_FLOW = true;
            mat->SW2_FLOW = true;
            mat->ivars.h   = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->ivars.uda = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->ivars.vda = mat->ivars.n; mat->ivars.n++; // 2D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '3') {
            mat->SW_FLOW = true;
            mat->SW3_FLOW = true;
            mat->ivars.dpl = mat->ivars.n; mat->ivars.n++; // 3D SW 
            mat->ivars.u   = mat->ivars.n; mat->ivars.n++; // 3D SW
            mat->ivars.v   = mat->ivars.n; mat->ivars.n++; // 3D SW
            mat->nSubmodels++;
        } else if (mat->code[0] == '4') {
            mat->NS_FLOW = true;
            mat->NS3_FLOW = true;
            mat->ivars.u   = mat->ivars.n;   mat->ivars.n++; // 3D NS 
            mat->ivars.v   = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->ivars.w   = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->ivars.prs = mat->ivars.n;   mat->ivars.n++; // 3D NS
            mat->nSubmodels++;
        } else if (mat->code[0] == '5') {
            mat->NS_FLOW = true;
            mat->NS3_SPLIT = true;
            mat->ivars.u = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->ivars.v = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->ivars.w = mat->ivars.n; mat->ivars.n++; // 3D NS SPLIT
            mat->nSubmodels++;
        } else if (mat->code[0] == '6') {
            mat->DW_FLOW = true;
            mat->ivars.h = mat->ivars.n; mat->ivars.n++;   // 2D OF
            mat->nSubmodels++;
        } else if (mat->code[0] == '7') {
            mat->WVEL_SPLIT = true;
            mat->ivars.w = mat->ivars.n; mat->ivars.n++;   // 3D SW - SPLIT
            mat->nSubmodels++;
        } else if (mat->code[0] == '8') {
            mat->PRESSURE = true;
            mat->ivars.prs = mat->ivars.n; mat->ivars.n++; // 3D NS - SPLIT
            mat->nSubmodels++;
        }

        // ground water
        if (mat->code[1] == '1') {
            mat->GW_FLOW = true;
            // only add depth if surface water has not
            if (mat->ivars.h == UNSET_INT) {
                mat->ivars.h = mat->ivars.n; mat->ivars.n++; // 1D/2D/3D GROUNDWATER
            }
            mat->nSubmodels++;
        }

        // transport
        mat->ntrns = 0;
        for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
            mat->TRANSPORT[itrns] = true;
            mat->ntrns++;
            mat->ivars.con[itrns] = mat->ivars.n; mat->ivars.n++; // 1D/2D/3D TRANSPORT
            mat->nSubmodels++;
        }
    assert(mat->nSubmodels > 0);
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
    if (mat->GW_FLOW)    {nSubMod_nvar[k]=1; k++;}
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
        if (mat->TRANSPORT[itrns]) {nSubMod_nvar[k]=1; k++;}    
    }

    int tmp[mat->nSubmodels], tmp2[mat->nSubmodels]; // Just to allocate
    sarray_init_value_int(tmp,mat->nSubmodels,UNSET_INT);
    sarray_init_value_int(tmp2,mat->nSubmodels,UNSET_INT);
    smodel_alloc_init_array(&(mat->model), tmp, tmp2, mat->nSubmodels, nSubMod_nvar); 

    // fill in model info
    int isubModel = 0;
    if (mat->SW1_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.uda;
        mat->model[isubModel].physics = SW1D_; // for body residuals
        mat->model[isubModel].physics_init = SW1D_; 
        isubModel++;
    }
    if (mat->SW2_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.uda;
        mat->model[isubModel].physics_vars[2] = mat->ivars.vda;
        mat->model[isubModel].physics = SW2D_;
        mat->model[isubModel].physics_init = SW2D_;
        isubModel++;
    }
    if (mat->SW3_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics_vars[1] = mat->ivars.u;
        mat->model[isubModel].physics_vars[2] = mat->ivars.v;
        mat->model[isubModel].physics = SW3D_;
        mat->model[isubModel].physics_init = SW3D_;
        isubModel++;
    }
    if (mat->NS3_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.u;
        mat->model[isubModel].physics_vars[1] = mat->ivars.v;
        mat->model[isubModel].physics_vars[2] = mat->ivars.w;
        mat->model[isubModel].physics_vars[3] = mat->ivars.prs;
        mat->model[isubModel].physics = NS3D_;
        mat->model[isubModel].physics_init = NS3D_;
        isubModel++;
    }
    if (mat->NS3_SPLIT) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.u;
        mat->model[isubModel].physics_vars[1] = mat->ivars.v;
        mat->model[isubModel].physics_vars[2] = mat->ivars.w;
        mat->model[isubModel].physics = NS3DSPLIT_;
        mat->model[isubModel].physics_init = NS3DSPLIT_;
        isubModel++;
    }
    if (mat->DW_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics = DW2D_;
        mat->model[isubModel].physics_init = DW2D_;
        isubModel++;
    }
    if (mat->WVEL_SPLIT) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.w;
        mat->model[isubModel].physics = WVEL_;
        mat->model[isubModel].physics_init = WVEL_;
        isubModel++;
    }
    if (mat->PRESSURE) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.prs;
        mat->model[isubModel].physics = PRS_;
        mat->model[isubModel].physics_init = PRS_;
        isubModel++;
    }
    if (mat->GW_FLOW) {
        mat->model[isubModel].physics_vars[0] = mat->ivars.h;
        mat->model[isubModel].physics = GW3D_;
        mat->model[isubModel].physics_init = GW3D_;
        isubModel++;
    }
    for (itrns=0; itrns<mat->ivars.ntrns; itrns++) {
        if (mat->TRANSPORT[itrns]) {
            mat->model[isubModel].physics_vars[0] = mat->ivars.con[itrns];
            mat->model[isubModel].physics = TRNS_;
            mat->model[isubModel].physics_init = TRNS_;
            isubModel++;
        }  
    }
    assert(isubModel == mat->nSubmodels);
}


