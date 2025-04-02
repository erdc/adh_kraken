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
int comm_update_ghost(double *arr_like,
                       SMPI *smpi
#ifdef _MESSG
                       ,MPI_Datatype data_type 
#endif
  )
{
#ifdef _MESSG
  //choose the neighborhood communicator
  //should we copy or just copy adress??
  MPI_Comm ADH_NEIGH = smpi->ADH_NEIGH;
  //for proof of concept, just nodal quantity for now

  //reallocate buffer
  messg_buffer_alloc(smpi->nrecv_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(data_type), /* the size of the items */
                    &(smpi->buffer_recv_neigh)   /* the buffer to be packed */
  );

  messg_buffer_alloc(smpi->nsend_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(data_type), /* the size of the items */
                    &(smpi->buffer_send_neigh)   /* the buffer to be packed */
  );


  //fill out recieving/sending buffer
  smpi->buffer_recv_neigh.buffer = &(arr_like[smpi->recv_ind]);
  double *bpntr  = (double *) smpi->buffer_send_neigh.buffer;


  for (int i = 0; i < smpi->nsend_neigh; i++ ){
    bpntr[i] = arr_like[smpi->dest_indices[i]];
  }


  int ierr = MPI_Neighbor_alltoallv(smpi->buffer_send_neigh.buffer, 
    smpi->dest_weights, smpi->dest_displs, data_type,
    smpi->buffer_recv_neigh.buffer, smpi->source_weights,
    smpi->source_displs, data_type, smpi->ADH_NEIGH);


  
#endif
  return ierr;
}
