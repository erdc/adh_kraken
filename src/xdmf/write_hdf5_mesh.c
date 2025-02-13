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

void write_hdf5_mesh(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt){
    //given a time step number, rights a mesh to the hdf5 file

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[RANK];
    hsize_t   offset[RANK];
    herr_t    status; 
    hsize_t   count[RANK];
    int nodes_per_pe,indx;
    int elem_per_pe=1;


    char number[50] = "";
    sprintf(number, "%d", nt);

     //mesh properties
    int numQuad;
    int numTri;
    int *connectivity;
    numQuad = 2;
    numTri = 2;




    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);

    //open file
    strcpy(fname,fbase);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);


    //group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Mesh/XY", H5P_DEFAULT);

    //////////////////////////////
    //Writing Nodes//////////////

    //declare global sizes of data set
    //Nodes first
    //this is local data
    nodes_per_pe=2;
    //xyz = (float *)malloc(sizeof(float) * nodes_per_pe * 2);
    //xyz = (float *)malloc(sizeof(float[nodes_per_pe][2]));
    float (*xyz)[2] = malloc(sizeof(float[nodes_per_pe][2]));

    //Hard code partition this time
    // 2 nodes per partition
    if(mpi_rank==0){
        xyz[0][0] = 0.0;xyz[0][1] = 0.0;
        xyz[1][0] = 1.0;xyz[1][1] = 0.0; }
    else if(mpi_rank==1){   
        xyz[0][0] = 1.0;xyz[0][1] = 1.0;
        xyz[1][0] = 0.0;xyz[1][1] = 1.0;} 
    else if(mpi_rank==2){   
        xyz[0][0] = 2.0;xyz[0][1] = 0.0;
        xyz[1][0] = 2.0;xyz[1][1] = 1.0; }
    else if(mpi_rank==3){
        xyz[0][0] = 2.0;xyz[0][1] = 2.0+nt; 
        xyz[1][0] = 0.0;xyz[1][1] = 2.0+nt;} 



    


    //create a dataspace
    // this is a global quantity
    dims[0]  = NumNodes;
    dims[1]  = 2;
    filespace = H5Screate_simple(2, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a simple dataspace where each process has a few of the rows
    count[0]  = nodes_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;

    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);


    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);


    //create data pointer
    //data = (int *)malloc(sizeof(int) * count[0] * count[1]);
    //for (i = 0; i < count[0] * count[1]; i++) {
    //    data[i] = mpi_rank + 10;
    //}

    //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xyz);
    free(xyz);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);
    //////////////////////////
    ////Nodal write complete//



    //Now write elemental data
    //////////////////////////
    //Element data set///////
    // Connectivity data 
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    grp2 = H5Gopen(file_id, "/Mesh/Elements", H5P_DEFAULT);




    //global attributes
    int nentry;
    nentry = numQuad*(4+1) + numTri*(3+1);
    
    //local data, different on each PE
    //first two PEs hold one quad. second PEs hold one triangle
    if(mpi_rank == 0 || mpi_rank ==1){ 
        connectivity = (int *)malloc(sizeof(int) * 5);}
    else if (mpi_rank ==2 || mpi_rank ==3){
        connectivity = (int *)malloc(sizeof(int) * 4);
    }
    
    //load in connectivity values
    if (mpi_rank==0){
        count[0] = 5;
        offset[0]=0;
        connectivity[0] = 5;
        connectivity[1] = 0;connectivity[2] = 1;connectivity[3] = 2;connectivity[4] = 3;   
    }
    else if(mpi_rank==1){
        count[0]=5;
        offset[0]=5;
        connectivity[0] = 5;
        connectivity[1] = 1;connectivity[2] = 4;connectivity[3] = 5;connectivity[4] = 2;
    }
    else if(mpi_rank==2){
        count[0]=4;
        offset[0]=10;
        connectivity[0] = 4;
        connectivity[1] = 2;connectivity[2] = 5;connectivity[3] = 6;   
        
    }
    else if(mpi_rank==3){
        count[0]=4;
        offset[0]=14;
        connectivity[0] = 4;
        connectivity[1] = 3;connectivity[2] = 2;connectivity[3] = 7;
    }

    //create dataspace and add
    //give global size
    dims[0] = nentry;
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp2, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, connectivity);
    free(connectivity);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp2);

    ////////////////////////////////////
    ///Elemental connections complete///

    //close mesh file
    H5Fclose(file_id);

}

