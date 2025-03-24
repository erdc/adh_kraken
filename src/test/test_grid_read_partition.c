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
    sgrid_read(&(dm.grid),dm.filename_grid,dm.grid->smpi->ADH_COMM);
#else
    sgrid_read(&(dm.grid),dm.filename_grid);
#endif

    sgrid_printScreen(dm.grid);

	return ierr;

}