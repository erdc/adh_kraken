/*! \file test_grid_read_partition.c This file tests the grid reader and partitioner*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests a basic sw2 case where we have sloping beach and still
 *  conditions. Mass conservation is tested
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int test_grid_read_partition(int npx, int npy, int nt){

	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// return error code
    //++++++++++++++++++++++++++++++++++++++++++++++
	int ierr = -1;
	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Create a design model
    //++++++++++++++++++++++++++++++++++++++++++++++
	SMODEL_DESIGN dm; smodel_design_defaults(&dm);
	
	//++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
	// Create a regular grid, from file
    //++++++++++++++++++++++++++++++++++++++++++++++
	dm.grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
	char filename[30];
    strcpy(filename,"../qa/flume/test"); printf("filename: %s\n",filename);
    strcpy(dm.filebase,filename);
    strcpy(filename,"../qa/flume/test.geo");
    strcpy(dm.filename_grid, filename);
    printf("filename: %s\n", dm.filename_grid);
#ifdef _MESSG
    MPI_Comm ADH_COMM;
    ADH_COMM = MPI_COMM_WORLD;
    sgrid_read(&(dm.grid),dm.filename_grid,ADH_COMM);
    ierr=0;
#else
    sgrid_read(&(dm.grid),dm.filename_grid);
#endif

    // all seems to work
    // in actual init
    // all other grid based fields will need to be filled here
    // then do a bandwidth minimizing partition
    ierr = comm_create_neighborhood(dm.grid);


    //see how to use neighborhood communicator
    double *arr = tl_alloc(sizeof(double), dm.grid->nnodes);
    sarray_init_value_dbl(arr, dm.grid->nnodes, dm.grid->smpi->myid);
#ifdef _MESSG
    ierr += comm_update_ghost(arr, dm.grid->smpi, MPI_DOUBLE);
#endif
    //after ghost update, all vals should match with .resident_pe
    //verify this condition
    for (int i = 0 ; i<dm.grid->nnodes ; i++){
        if (arr[i] != dm.grid->node[i].resident_pe){ierr+=1;}
    }


    int *iarr = tl_alloc(sizeof(double), dm.grid->nnodes);
    sarray_init_value_int(iarr, dm.grid->nnodes, dm.grid->smpi->myid);
#ifdef _MESSG
    ierr += comm_update_ghost(iarr, dm.grid->smpi, MPI_INT);
#endif
    for (int i = 0 ; i<dm.grid->nnodes ; i++){
        if (iarr[i] != dm.grid->node[i].resident_pe){ierr+=1;}
    }

    //sgrid_printScreen(dm.grid);

    //now to comply with unit tests either return 1 if fail or 0 if no fail
    if(ierr!=0){ierr=1;}

	return ierr;

}