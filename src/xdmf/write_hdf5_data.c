#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>
#define NumElements 4
#define NumNodes 8
#define NEntry 18
#define RANK 2
//RANK is dimension of data set, 2 is a matrix, 1 is a vector

void write_hdf5_data(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt, float t)
{
    //writes hdf5 data at particular time step nt
    //for this simple test case assert the size is only 4
    assert(mpi_size==4);


    hid_t     file_id, dset_id, memspace, filespace, plist_id, grp;
    hsize_t   dims[RANK];
    //float *xyz;
    int i;
    hsize_t   offset[RANK];
    herr_t    status; 
    hsize_t   count[RANK];
    int nodes_per_pe,indx;
    nodes_per_pe=2;
    int elem_per_pe=1;


    //mesh properties
    int numQuad;
    int numTri;
    int *connectivity;
    numQuad = 2;
    numTri = 2;


    char number[50] = "";
    sprintf(number, "%d", nt);

    //tutorial at
    //https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html



    //Open file for parallel i/o
    ////////////////////////////////////////////////
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //open file
    char fname[50];
    strcpy(fname,fbase);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    ////////////////////////////////////////////////



    
    //Open up data groups
    //group 2 contains any info with data
    grp = H5Gopen(file_id, "/Data/NodalScalar", H5P_DEFAULT);




    
    //////Scalar data write////////////////
    //Nodally based data first

    float *scalardata; 
    scalardata = (float *)malloc(sizeof(float) * nodes_per_pe);
        
    for(indx = 0; indx < nodes_per_pe; indx++)
    {   
        scalardata[indx] = (mpi_rank*nodes_per_pe+indx)*100 + t*50;
    }

    dims[0] = NumNodes;
    count[0] = nodes_per_pe;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(1, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, scalardata);
    free(scalardata);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);



    //End of Nodal based scalar

    


    //Write some elementally based data ////////////
    //Elementally based data
    grp = H5Gopen(file_id, "/Data/ElementalScalar", H5P_DEFAULT);
     
    float *scalardata2; 
    scalardata2 = (float *)malloc(sizeof(float) * elem_per_pe);
        
    for(indx = 0; indx < 1; indx++)
    {   
        scalardata2[indx] = (mpi_rank)*100 + t*50;
    }

    dims[0] = NumElements;
    count[0] = 1;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(1, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, scalardata2);
    free(scalardata2);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);

    //End of elementally based data
    /////////////////////////////////////////////

    
    ////////////////////////////////////////////
    
    //vector based nodal data
    grp = H5Gopen(file_id, "/Data/NodalVector", H5P_DEFAULT);
    float (*vectordata)[2] = malloc(sizeof(float[nodes_per_pe][2]));
        
    for(indx = 0; indx < nodes_per_pe; indx++)
    {   
        vectordata[indx][0] = 1+t;
        vectordata[indx][1] = 1;
    }

    dims[0] = NumNodes;
    dims[1] = 2;
    count[0]  = nodes_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    filespace = H5Screate_simple(2, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, vectordata);
    free(vectordata);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);
    //end of nodal vectors
    /////////////////////////////////////////////////////


    //Elementally based vectors
    //////////////////////////////////////////////////////
    grp =  H5Gopen(file_id, "/Data/ElementalVector", H5P_DEFAULT);
    float (*vectordata2)[2] = malloc(sizeof(float[elem_per_pe][2]));
        
    for(indx = 0; indx < elem_per_pe; indx++)
    {   
        vectordata2[indx][0] = 1-t;
        vectordata2[indx][1] = 1;
    }

    dims[0] = NumElements;
    dims[1] = 2;
    count[0]  = elem_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    filespace = H5Screate_simple(2, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, vectordata2);
    free(vectordata2);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);

    



    /////////////////////////////////////////////////////
    H5Fclose(file_id);

    if (mpi_rank == 0)
        printf("PHDF5 example finished with no errors\n");
  
}
