#ifndef _H_MODELS_
#define _H_MODELS_
//Integer codes for routine pointers


/* Variable codes */
//Dont think we will need these any longer
//#define PERTURB_U 1
//#define PERTURB_V 2
//#define PERTURB_W 3
//#define PERTURB_DPL 4
//#define PERTURB_H 5
//#define PERTURB_C 6
//#define PERTURB_P 7
//#define PERTURB_D 8
//#define PERTURB_NONE -1




//// Definitions for residual/init function pointer calls
//// to simplify code, must always be contiguous
//// that way we can use array look up as map
//#define N_RESID_ROUTINES 40
//#define NO_RESID_ 0
//#define SW1D_ 1
//#define SW2D_ 2
//#define SW3D_ 3
//#define WVEL_ 4
//#define NS1D_ 5
//#define NS2D_ 6
//#define NS3D_ 7
//#define NS1DSPLIT_ 8
//#define NS2DSPLIT_ 9
//#define NS3DSPLIT_ 10
//#define PRS_  11
//#define DW1D_ 12
//#define DW2D_ 13
//#define GW1D_ 14
//#define GW2D_ 15
//#define GW3D_ 16
//#define TRNS_ 17
//#define TRNS1D_ 18
//#define TRNS2D_ 19
//#define TRNS3D_ 20
//#define HEAT1D_ 21
//#define HEAT2D_ 22
//#define HEAT3D_ 23
//#define POISSON2D_ 24
//// Definitions for boundary treatments, included as resid routines
//#define SW2_BC_DISCHARGE_ 25
//#define SW2_BC_ELE_ 26
//#define SW2_BC_FLAPD_ 27
//#define SW2_BC_FLAPU_ 28
//#define SW2_BC_FLUX_ 29
//#define SW2_BC_H_ 30
//#define SW2_BC_HYBRID_ 31
//#define SW2_BC_OUTFLOW_ 32
//#define SW2_BC_SLUICED_ 33
//#define SW2_BC_SLUICEU_ 34
//#define SW2_BC_TAILWATER_ 35
//#define SW2_BC_VEL_ 36
//#define SW2_BC_VEL_ELE_ 37
//#define SW2_BC_WEIRD_ 38
//#define SW2_BC_WEIRU_ 39//
//

////init routine codes
//#define N_INIT_ROUTINES 2
//#define NO_INIT 0
//#define SW2_INIT 1

//time step routines
//forward step pointer codes
//see time_loop
#define N_TIME_STEPPERS 1
#define FE_NEWTON 0



// Definitions for residual/init function pointer calls
//#define N_RESID_ROUTINES 25
//#define NO_RESID_ 0
//#define SW1D_ 1
//#define SW2D_ 2
//#define SW3D_ 3
//#define WVEL_ 31
//#define NS1D_ 11
//#define NS2D_ 12
//#define NS3D_ 13
//#define NS1DSPLIT_ 21
//#define NS2DSPLIT_ 22
//#define NS3DSPLIT_ 23
//#define PRS_  4
//#define DW1D_ 51
//#define DW2D_ 52
//#define GW1D_ 61
//#define GW2D_ 62
//#define GW3D_ 63
//#define TRNS_ 70
//#define TRNS1D_ 71
//#define TRNS2D_ 72
//#define TRNS3D_ 73
//#define HEAT1D_ 81
//#define HEAT2D_ 82
//#define HEAT3D_ 83
//#define POISSON2D_ 92
//// Definitions for boundary treatments, included as resid routines
//#define SW2_BC_DISCHARGE_ 200
//#define SW2_BC_ELE_ 201
//#define SW2_BC_FLAPD_ 202
//#define SW2_BC_FLAPU_ 203
//#define SW2_BC_FLUX_ 204
//#define SW2_BC_H_ 205
//#define SW2_BC_HYBRID_ 206
//#define SW2_BC_OUTFLOW_ 207
//#define SW2_BC_SLUICED_ 208
//#define SW2_BC_SLUICEU_ 209
//#define SW2_BC_TAILWATER_ 210
//#define SW2_BC_VEL_ 211
//#define SW2_BC_VEL_ELE_ 212
//#define SW2_BC_WEIRD_ 213
//#define SW2_BC_WEIRU_ 214
////init routine codes
//#define N_INIT_ROUTINES 2
//#define NO_INIT 0
//#define SW2 1





#endif
