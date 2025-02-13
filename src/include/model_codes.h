//Integer codes for routine pointers


/* Variable codes */
#define PERTURB_U 1
#define PERTURB_V 2
#define PERTURB_W 3
#define PERTURB_DPL 4
#define PERTURB_H 5
#define PERTURB_C 6
#define PERTURB_P 7
#define PERTURB_D 8
#define PERTURB_NONE -1


//forward step pointer codes
//see time_loop
#define N_TIME_STEPPERS 1
#define FE_NEWTON 0

//residual routine codes
#define N_RESID_ROUTINES 3
#define SW2 0
#define POISSON 1
#define HEAT 2

//init routine codes
#define N_INIT_ROUTINES 1
#define SW2 0

// Definitions for residual/init function pointer calls
#define SW1D_ 1
#define SW2D_ 2
#define SW3D_ 3
#define WVEL_ 31
#define NS1D_ 11
#define NS2D_ 12
#define NS3D_ 13
#define NS1DSPLIT_ 21
#define NS2DSPLIT_ 22
#define NS3DSPLIT_ 23
#define PRS_  4
#define DW1D_ 51
#define DW2D_ 52
#define GW1D_ 61
#define GW2D_ 62
#define GW3D_ 63
#define TRNS_ 70
#define TRNS1D_ 71
#define TRNS2D_ 72
#define TRNS3D_ 73
#define HEAT1D_ 81
#define HEAT2D_ 82
#define HEAT3D_ 83


//global array of function pointers to resid routines
int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int);
//global array of function pointers to init routines
int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *);


