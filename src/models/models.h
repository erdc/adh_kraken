#ifndef H_MODELS_
#define H_MODELS_


void set_function_pointers( int (*resid_routines[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int), int (*init_routines[N_INIT_ROUTINES])(SMODEL_SUPER *));


#endif
