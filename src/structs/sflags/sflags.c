/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sflags.c This file collects methods of the SFLAGS structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Intializes a SFLAGS structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] f (SFLAGS *)  a pointer to an AdH superModel SFLAGS structure
 *
 * \note CJT: Some of these can go
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sflags_init(SFLAGS *f) {
    
    f->NS_FLOW = false;
    f->NS2_FLOW = false;
    f->NS3_FLOW = false;
    f->GW_FLOW = false;
    f->SW_FLOW = false;
    f->SW2_FLOW = false;
    f->SW3_FLOW = false;
    f->DIFFUSIVE_WAVE = false;
    f->MG = false;
    f->BAROCLINIC = false;
    f->EOS = 0;
    f->GRID_ADAPTION = false;
    f->ADAPTED_THE_GRID = false;
    f->GRID_REFINED = false;
    f->GRID_UNREFINED = false;
    f->TIME_ADAPT = false;
    f->TIME_ADAPT_FAIL = false;
    f->TRANSPORT = false;
    f->NS2_TRANSPORT = false;
    f->NS3_TRANSPORT = false;
    f->SW2_TRANSPORT = false;
    f->SW3_TRANSPORT = false;
    f->VORTICITY = false;
    f->SEDIMENT = false;
    f->SEDLIB = false;
#ifdef _ADH_GROUNDWATER
    f->GW_SALINITY = false;
    f->GW_TRANSPORT = false;
    f->GW_REACTION = false;
    f->RAY_TRACING = false;
    f->SOCKETS = false;
    f->METEOROLOGY = false;
#endif
    f->STEADY_STATE = false;
    f->ICM = false;
    f->NSM = false;
    f->WAVE_RADS = false;
    f->WAVE_STRESS = false;
    f->WIND_STRESS = false;
    f->WAVE_STATION = false;
    f->WIND_STATION = false;
    f->CORIOLIS = false;
    f->MUC = false;
    f->UNITS = false;
    f->CONVEYANCE = false;
    f->PRN_ADPT = false;
    f->OUTPUT = false;
    f->SOLVE_ATF = false;
    f->UMFPACK_FAIL = false;
    f->ICE = false;
    f->INS = false;
    f->TIDE = false;
    f->UNREFINE = false;
    f->CSTORM_WSID = false;
    f->PC_FILE_XDMF = false;
    f->FLUX_WEIGHTED_NORMALS = false;
    
}
