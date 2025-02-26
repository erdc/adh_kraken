#ifndef H_SMAT_PHYSICS_
#define H_SMAT_PHYSICS_

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    int id;
    char code[10];
    int ntrns; // # of transport constituents on this element
    SIVAR_POSITION ivar_pos; // a list of independent variables on the element
    int n; //number of variables
    int *ivar_loc; //arry of n who's entries are location in **ivar array
    int nSubmodels;
    SMODEL *model; // [nSubModels] length array

    bool SW_FLOW;   // 1, 2, 3
    bool SW1_FLOW;  // 1
    bool SW2_FLOW;  // 2
    bool SW3_FLOW;  // 3
    bool NS_FLOW;   // 4
    bool NS3_FLOW;  // 4
    bool NS3_SPLIT; // 5
    bool DW_FLOW;   // 6
    bool WVEL_SPLIT;// 7
    bool PRESSURE;  // 8
    bool GW_FLOW;   // 9
    bool VORTICITY;
    bool SEDIMENT;
    bool SEDLIB;
    bool ICM;
    bool NSM;
    bool WAVE;
    bool WIND;
    bool *TRANSPORT; // ntrans

} SMAT_PHYSICS;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// Methods
void smat_physics_alloc_init(SMAT_PHYSICS *mat, char codes[10]);
void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, char codes[][10]);
void smat_physics_alloc_init_ptr(SMAT_PHYSICS *mat, char *codes);
void smat_physics_alloc_init_ptr_array(SMAT_PHYSICS **mat_physics, int nmat, char **codes);
void smat_physics_free_array(SMAT_PHYSICS *mat, int nmat);
void smat_physics_position_flag(SMAT_PHYSICS **mat_node, int nnodes, int *FLAG);
void smat_physics_printScreen(SMAT_PHYSICS *m);
void smat_physics_update_array(SMAT_PHYSICS *m, int nmat_physics, SIVAR_POSITION* ivar_pos);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
