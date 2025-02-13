#ifndef H_SCOVERAGE_
#define H_SCOVERAGE_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
   \file type_scoverage.h
   \brief Data Structure for AdH Physics/Parameter coverages
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {

   int nelems1D,nelems2D,nelems3D; // for convienence 
   int ncoverages1D;
   int ncoverages2D;
   int ncoverages3D;
   int **coverage_1D;
   int **coverage_2D;
   int **coverage_3D;

} SCOVERAGE;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Struct Methods
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void scoverage_alloc_init(SCOVERAGE **coverage,int nelems1D,int nelems2D,int nelems3D,int ncoverages1D,int ncoverages2D,int ncoverages3D);
void scoverage_free(SCOVERAGE *c);
void scoverage_read(SCOVERAGE **c,char *filebase);


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#endif