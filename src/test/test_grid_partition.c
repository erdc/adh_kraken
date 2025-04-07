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
int test_grid_partition(int comm_type){

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
    ierr = comm_create_neighborhood(dm.grid, comm_type);


    ierr = sgrid_repartition(dm.grid);
    //sgrid_printScreen(dm.grid);
    //now to comply with unit tests either return 1 if fail or 0 if no fail
    if(ierr!=0){ierr=1;}

    //free stuff, need to fix for parallel
    smodel_design_free(&dm);

	return ierr;

}