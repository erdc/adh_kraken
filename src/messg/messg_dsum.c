/*! \file  messg_dsum.c This file has function that computes sum of value over all processors*/
#include "adh.h"
/*!
   \brief Sums a single double value across processors

   \param[in,out] x (double) - Value in
   \returns Sum of Value Across Processors
 */
double messg_dsum(double x
  #ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
                  )
{
#ifdef _MESSG
  double x_send = DZERO;        /* the variables for the pass */
  int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

  x_send = x;
  ierr_code = MPI_Allreduce(&x_send, &x, 1, MPI_DOUBLE, MPI_SUM, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
#endif
  return (x);
}
