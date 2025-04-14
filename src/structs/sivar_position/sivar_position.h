#ifndef H_SIVAR_POS_
#define H_SIVAR_POS_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 \file sivar_pos.h
 \brief Stores the independent variable positions in the **var storage
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    
    int n;     // total # of independent variables active (includes transport)
    int ntrns; // total # of transport variables active only
    int *var;  // position in array form of length adh_def.n_ivars ( intended
               // to be used as map that takes in index of adh_def and produces index number
               // in the relevant structure (mat_physics or super model))
    int *var_code;  //array of ivar_pos.n who's entries are location of independent variables
                    // in adh_defs (Not supermodel), used in assemble_jacobian
                    // should/do we need it in SIVAR_POSITION (remember there is a member in 
                    // smat physics and smodel_super, in smodel_super this array would tell
                    // us the location of indpendent variables in adh_defs for all active variabses
                    // in the supermodel)
                    //NOTE: can easilt construct this by looping through ivar_pos.var and
                    //putting any non UNSET_INTs in contiguous entries
                    //This is the reverse of the map in *var, takes in local variable number
                    // (whether in physics mat or super model) and gives variable number
                    // in adh_def
    
} SIVAR_POSITION;   

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void sivar_position_init(SIVAR_POSITION *ip);
void sivar_position_free(SIVAR_POSITION *ip);
void sivar_position_printScreen(SIVAR_POSITION *ip);
void sivar_position_map(SIVAR_POSITION *ip, int *FLAGS);
int  sivar_position_build_dof_map(SIVAR_POSITION *ip, int nnodes,  SIVAR_POSITION *ipNode, int ***ivars);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif





