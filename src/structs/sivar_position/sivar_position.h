#ifndef H_SIVAR_POS_
#define H_SIVAR_POS_

#define MAX_VARS 10 // add to this if you add independent variables
#define MAX_TRNS_VARS 20 // maximum # of transport independent variables. Can be increased.

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

   // independent variable order!
   int h;
   int u;
   int v;
   int w;
   int uda;
   int vda;
   int dpl;
   int prs;
   int heat;
   int sal;
   int con[MAX_TRNS_VARS];

} SIVAR_POSITION;   

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void sivar_position_init(SIVAR_POSITION *ip);
void sivar_position_printScreen(SIVAR_POSITION *ip);
void sivar_position_map(SIVAR_POSITION *ip, int *FLAGS);
int sivar_position_build_dof_map(SIVAR_POSITION *ip, int nnodes,  SIVAR_POSITION *ipNode, int ***ivars);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif



