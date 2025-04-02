#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Updates ghost dofs by sending values from owning process
 *             for double, using neighborhood AlltoAllv
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 * \note The general routine should update ghost values on each process from the owners
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//void comm_update_ghost(void *arr_like,
//                       SMPI *smpi,
//                       SGRID *grid
//#ifdef _MESSG
//                       , MPI_Datatype data_type 
//#endif
//  )
void comm_update_double(double *vec, int size_v, int npe, int rank)
{
//#ifdef _MESSG
//  //choose the neighborhood communicator
//  //should we copy or just copy adress??
//  MPI_Comm ADH_NEIGH = smpi->ADH_NEIGH;
//  //for proof of concept, just nodal quantity for now
//  int nghost = grid->nnodes - grid->my_nnodes;//
//
//

//  //for now realloc buffer stuff in here, need to put somewhere else
//  messg_buffer_alloc(nghost,  /* the number of items to be stored in the buffer */
//                    sizeof(data_type), /* the size of the items */
//                    smpi->buffer_recv_neigh   /* the buffer to be packed */
//  );
//  //simple double for testing then will be the weights from prior//
//
//

//  messg_buffer_alloc(nghost,  /* the number of items to be stored in the buffer */
//                    sizeof(data_type), /* the size of the items */
//                    smpi->buffer_send_neigh   /* the buffer to be packed */
//  );//
//

//  //for now make simple reciever
//  int *rdispls = (int*) tl_alloc(sizeof(int),smpi->indegree);
//  int ctr=0;
//  for(int i = 0; i < smpi->indegree; i++ ){
//    rdispls[i] = ctr;
//    ctr+=grid->smpi->source_weights[i];
//  }//

//  double *recvbuf = (double *) tl_alloc(sizeof(double), ctr);//
//
//

//  int ierr = MPI_Neighbor_alltoallv(arr_like, const int sendcounts[],
//    const int sdispls[], data_type,
//    recvbuf, grid->smpi->source_weights,
//    rdispl, data_type, smpi->ADH_NEIGH);
//  
//#endif
  return;
}
