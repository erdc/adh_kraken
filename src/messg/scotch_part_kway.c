#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Redistributes nodes to minimize bandwidth using PTSCOTCH library
 * 			   Small wrapper that fills the *part arry
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] pgrid (SGRID *)  pointer to an AdH grid
 * @returns int - 0 if succesful
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int scotch_partkway(int local_nnode, int *xadj, int *adjncy, int *vwgt, int *adjwgt, 
	int stratflag, int numflag, int npes, int balance_ratio, int *part, MPI_Comm *comm){
	int ierr = 0;
#ifdef _MESSG
#ifdef _SCOTCH
	//quick check that int aligns with SCOTCH's wrapper
	assert(sizeof(SCOTCH_Num) == sizeof(int));
	assert(SCOTCH_numSizeof () == sizeof (SCOTCH_Num));

	//create a SCOTCH distributed graph
	SCOTCH_Dgraph grafdat;

	ierr = SCOTCH_dgraphInit(&grafdat, *comm);
	if (ierr!=0){return -1;}

	//build the distributed graph
	SCOTCH_Num baseval = numflag;
	assert(baseval == 0); //should always be 0 based
	SCOTCH_Num vertlocnbr = local_nnode;
	SCOTCH_Num edgelocnbr = xadj[vertlocnbr] - baseval;
	SCOTCH_Strat stradat;
	ierr = SCOTCH_dgraphBuild(&grafdat, baseval, vertlocnbr, vertlocnbr,
		xadj, xadj+1, vwgt, NULL, edgelocnbr, edgelocnbr, adjncy, NULL, adjwgt);
	if (ierr!=0){return ierr;}
		
	ierr = SCOTCH_stratInit(&stradat);
	if (ierr!=0){return ierr;}
	ierr = SCOTCH_stratDgraphMapBuild(&stradat, stratflag, npes, npes, balance_ratio);
	if (ierr!=0){return ierr;}
	
	//should we call DGraphMap or dgraphPart??
	//Idk which one is better
	ierr = SCOTCH_dgraphPart(&grafdat, npes, &stradat, part);
	//alternative using MapBuild
	//SCOTCH_Arch archdat;
	//ierr = SCOTCH_archInit(&archdat);
	//ierr = SCOTCH_archCmplt(&archdat, npes);
	//ierr = SCOTCH_dgraphMap(&grafdat, &archdat &stradat, part);
	//ierr = SCOTCH_archExit(&archdat);
	
	SCOTCH_stratExit(&stradat);
	SCOTCH_dgraphExit(&grafdat);
#endif
#endif
	return  ierr;
}