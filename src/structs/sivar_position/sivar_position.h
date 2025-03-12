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
    int var[N_IVARS_TOTAL]; // position in array form
    
} SIVAR_POSITION;   

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void sivar_position_init(SIVAR_POSITION *ip);
void sivar_position_printScreen(SIVAR_POSITION *ip);
void sivar_position_map(SIVAR_POSITION *ip, int *FLAGS);
int  sivar_position_build_dof_map(SIVAR_POSITION *ip, int nnodes,  SIVAR_POSITION *ipNode, int ***ivars);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif





