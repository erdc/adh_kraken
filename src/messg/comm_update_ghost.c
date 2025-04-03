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
int comm_update_ghost(void *arr_like,
                       SMPI *smpi
#ifdef _MESSG
                       ,MPI_Datatype data_type 
#endif
  )
{
#ifdef _MESSG

  if(smpi->npes == 1){return 0;}

  //reallocate buffer
  if (data_type == MPI_DOUBLE){
    messg_buffer_alloc(smpi->nrecv_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(double), /* the size of the items */
                    &(smpi->buffer_recv_neigh)   /* the buffer to be packed */
    );
    messg_buffer_alloc(smpi->nsend_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(double), /* the size of the items */
                    &(smpi->buffer_send_neigh)   /* the buffer to be packed */
    );
  }else if(data_type == MPI_INT){
    messg_buffer_alloc(smpi->nrecv_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(int), /* the size of the items */
                    &(smpi->buffer_recv_neigh)   /* the buffer to be packed */
    );
    messg_buffer_alloc(smpi->nsend_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(int), /* the size of the items */
                    &(smpi->buffer_send_neigh)   /* the buffer to be packed */
    );

  }
  
  //typecast the vector of intrest based on data_type
  //typecast the send buffer and set values
  double *ptr_double;
  double *bpntr_double;
  int *ptr_int;
  int *bpntr_int;

  if (data_type == MPI_DOUBLE){
    ptr_double = (double *) arr_like;
    bpntr_double  = (double *) smpi->buffer_send_neigh.buffer;
    //recieving buffer, assumed contiguous and in order
    smpi->buffer_recv_neigh.buffer =  &ptr_double[smpi->recv_ind];
    for (int i = 0; i < smpi->nsend_neigh; i++ ){
      //fill out sending buffer
      bpntr_double[i] =  ptr_double[smpi->dest_indices[i]];
    }
  }else if (data_type == MPI_INT){
    ptr_int = (int *) arr_like;
    bpntr_int  = (int *) smpi->buffer_send_neigh.buffer;
    //recieving buffer, assumed contiguous and in order
    smpi->buffer_recv_neigh.buffer =  &ptr_int[smpi->recv_ind];
    for (int i = 0; i < smpi->nsend_neigh; i++ ){
      //fill out sending buffer
      bpntr_int[i] =  ptr_int[smpi->dest_indices[i]];
    }
  }
  

 

  int ierr = MPI_Neighbor_alltoallv(smpi->buffer_send_neigh.buffer, 
    smpi->dest_weights, smpi->dest_displs, data_type,
    smpi->buffer_recv_neigh.buffer, smpi->source_weights,
    smpi->source_displs, data_type, smpi->ADH_NEIGH);


  
#endif
  return ierr;
}
