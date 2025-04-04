#ifndef H_MESSG_
#define H_MESSG_

//MPI related routines
double messg_dmax(double x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
  );

double messg_dmin(double x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
                  );

int messg_imin(int x            /* the value */
#ifdef _MESSG
               , MPI_Comm ADH_COMM
#endif
  );

int messg_imax(int x            /* the value */
#ifdef _MESSG
               , MPI_Comm ADH_COMM
#endif
  );

double messg_dsum(double x
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
                  );

void messg_err(int ierr);

void messg_buffer_init(MESSG_BUFFER * buffer,   /* the message buffer */
                       int i_processor  /* the target processor */
  );

void messg_buffer_free(MESSG_BUFFER * buffer    /* the buffer to be packed */
  );

void messg_gather_int(int root,
                      int *my_x, /* my part of the array */
                      int *x,    /* the resulting array */
                      int size  /* the size of x */
#ifdef _MESSG
                      ,MPI_Comm ADH_COMM
#endif
    );
void messg_gather_dbl(int root,
                      double *my_x, /* my part of the array */
                      double *x,    /* the resulting array */
                      int size  /* the size of x */
#ifdef _MESSG
                      ,MPI_Comm ADH_COMM
#endif
    );


int comm_update_double(double *arr, SMPI *smpi, int type);
int comm_update_int(int *arr,SMPI *smpi, int type);

#ifdef _MESSG
void messg_barrier(MPI_Comm ADH_COMM);
#endif

#ifdef _MESSG
    int messg_comm_rank(MPI_Comm comm);
#else
    int messg_comm_rank(void);
#endif

#ifdef _MESSG
int messg_comm_size(MPI_Comm ADH_COMM);
#else
int messg_comm_size(void);
#endif

#ifdef _MESSG
void messg_arecv(MESSG_BUFFER *buffer, /* the message */
                 int tag_1,        /* the tag */
                 SMPI *smpi
  );
#endif

void messg_buffer_alloc(int nitem,  /* the number of items to be stored in the buffer */
                        size_t sizeof_item, /* the size of the items */
                        MESSG_BUFFER * buffer   /* the buffer to be packed */
  );

void messg_asend(MESSG_BUFFER * buffer, /* the message */
                 int tag,        /* the tag */
                 SMPI *smpi
  );

void messg_wait(SMPI *smpi);



int comm_create_neighborhood(SGRID *grid, int type);
#endif
