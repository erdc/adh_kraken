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
int comm_update_int(int *arr,SMPI *smpi, int type)
{
        int ierr = 0;
#ifdef _MESSG



  if(smpi->npes == 1){return 0;}

  messg_buffer_alloc(smpi->nrecv_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(int), /* the size of the items */
                    &(smpi->buffer_recv_neigh)   /* the buffer to be packed */
  );
  messg_buffer_alloc(smpi->nrecv_neigh,  /* the number of items to be stored in the buffer */
                    sizeof(int), /* the size of the items */
                    &(smpi->buffer_send_neigh)   /* the buffer to be packed */
  );


    //typecast the vector of intrest based on data_type
    //typecast the send buffer and set values
    int *ptr_int;
    int *bpntr_int;
    int k = 0;
    ptr_int = (int *) arr;
    bpntr_int  = (int *) smpi->buffer_send_neigh.buffer;

  if (type == NEIGHBOR){


    //recieving buffer, assumed contiguous and in order
    smpi->buffer_recv_neigh.buffer =  &ptr_int[smpi->recv_ind];
    for (int i = 0; i < smpi->nsend_neigh; i++ ){
      //fill out sending buffer
      bpntr_int[i] =  ptr_int[smpi->dest_indices[i]];
    }


    ierr = MPI_Neighbor_alltoallv(smpi->buffer_send_neigh.buffer, 
      smpi->dest_weights, smpi->dest_displs, MPI_INT,
      smpi->buffer_recv_neigh.buffer, smpi->source_weights,
      smpi->source_displs, MPI_INT, smpi->ADH_NEIGH);
  
  }else if (type == P2P){

    //recieving buffer, assumed contiguous and in order
    for(int i = 0 ; i <smpi->indegree ; i++){
        smpi->buffer_recv_neigh.buffer =  &ptr_int[smpi->recv_ind + smpi->source_displs[i]];
        ierr = MPI_Irecv(smpi->buffer_recv_neigh.buffer,
                  smpi->source_weights[i], MPI_INT, smpi->sources[i], TAG_UPDATE,
                  smpi->ADH_COMM, &(smpi->msg_request[i]));
    }


    for (int i = 0; i < smpi->outdegree; i++){

        for (int j = 0; j < smpi->dest_weights[i]; j++ ){
            //fill out sending buffer
            bpntr_int[j] =  ptr_int[smpi->dest_indices[k]];
            k++;
        }

        ierr = MPI_Isend(smpi->buffer_send_neigh.buffer,
                  smpi->dest_weights[i], MPI_INT, smpi->dest[i], TAG_UPDATE,
                  smpi->ADH_COMM, &(smpi->msg_request[i+smpi->indegree]));

    }

      ierr = MPI_Waitall(smpi->outdegree, smpi->msg_request, MPI_STATUS_IGNORE);


  }else{
    //return error code if option isnt correct
    ierr = 1;
  }

#endif
  return ierr;
}