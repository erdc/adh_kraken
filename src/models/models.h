#ifndef _H_MODELS_
#define _H_MODELS_



void set_models(int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int), int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *));
//global array of function pointers to resid routines
int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int);
//global array of function pointers to init routines
int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *);
#endif
