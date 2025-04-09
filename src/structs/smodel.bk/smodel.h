#ifndef H_SMODEL_
#define H_SMODEL_


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    int physics;        // a model ID for elemental residual functions, etc.
    int physics_init;   // a model ID for elemental body initialization
    int nvar;           // the # of independent variable on this model
    int *physics_vars;  // the location of the variables in the superModel 
    
} SMODEL;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void smodel_alloc_init(SMODEL *m, int phys, int physInit, int nvar);
void smodel_alloc_init_array(SMODEL **mod, int *phys, int *physInit, int nmods, int *nvars);
void smodel_free(SMODEL *model);
void smodel_free_array(SMODEL *mods, int nmods);
void smodel_printScreen(SMODEL *m);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
