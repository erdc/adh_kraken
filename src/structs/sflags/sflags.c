/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sflags.h This file collects methods of the SFLAGS structure */
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

  f->NS_FLOW = OFF;
  f->NS2_FLOW = OFF;
  f->NS3_FLOW = OFF;
  f->GW_FLOW = OFF;
  f->SW_FLOW = OFF;
  f->SW2_FLOW = OFF;
  f->SW3_FLOW = OFF;
  f->DIFFUSIVE_WAVE = OFF;
  f->MG = OFF;
  f->BAROCLINIC = OFF;
  f->EOS = 0; 
  f->GRID_ADAPTION = OFF;
  f->ADAPTED_THE_GRID = NO;
  f->GRID_REFINED = NO;
  f->GRID_UNREFINED = NO;
  f->TIME_ADAPT = OFF;
  f->TIME_ADAPT_FAIL = NO;
  f->TRANSPORT = OFF;
  f->NS2_TRANSPORT = OFF;
  f->NS3_TRANSPORT = OFF;
  f->SW2_TRANSPORT = OFF;
  f->SW3_TRANSPORT = OFF;
  f->VORTICITY = OFF;
  f->SEDIMENT = OFF;
  f->SEDLIB = OFF;
  // f->GW_SALINITY = OFF;
  // f->GW_TRANSPORT = OFF;
  // f->GW_REACTION = OFF;
  // f->RAY_TRACING = OFF;
  // f->SOCKETS = OFF;
  // f->METEOROLOGY = OFF;
  f->STEADY_STATE = OFF;
  f->ICM = OFF;
  f->NSM = OFF;
  f->WAVE = OFF;
  f->WIND = OFF;
  f->WAVE_STATION = OFF;
  f->WIND_STATION = OFF;
  f->WIND_LIBRARY = OFF;
  f->CORIOLIS = OFF;
  f->MUC = OFF;
  f->UNITS = OFF;
  f->CONVEYANCE = OFF;
  f->PRN_ADPT = OFF;
  f->OUTPUT = OFF;
  f->SOLVE_ATF = OFF;
  f->UMFPACK_FAIL = NO;
  f->ICE = OFF;
  f->INS = OFF;
  f->TIDE = OFF;
  f->UNREFINE = NO;
  f->CSTORM_WSID = OFF;
  f->PC_FILE_XDMF = OFF;
  f->FLUX_WEIGHTED_NORMALS = OFF;
  f->fout_bed_velocity = OFF;
  f->fout_surface_velocity = OFF;
  f->fout_depth_avg_velocity = OFF;
  f->fout_pressure = OFF;
  f->fout_grid_speed = OFF;
  f->fout_wind = OFF;
  f->fout_wave = OFF;
  f->fout_hyd_vis = OFF;
  f->fout_trn_dif = OFF;
  f->fout_chop = OFF;
  f->fout_grid2dm = OFF;
  f->fout_adaption = OFF;
  f->fout_adapt_grid = OFF;
  f->fout_adapt_sw = OFF;
  f->fout_adapt_ns = OFF;
  f->fout_adapt_con = OFF;
  f->fout_adapt_sed = OFF;
  f->fout_adapt_gw = OFF;
  
}