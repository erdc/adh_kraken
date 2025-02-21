#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
void hdf5_init(SGRID *g, char *fbase) {
    
#ifdef _HDF5
    
    //++++++++++++++++++++++++++++++++++++++
    // Open HDf5 File and Initiate
    //++++++++++++++++++++++++++++++++++++++
    // 1. Create HDf5 File for the simulation
    // 2. Create all necessary data groups for the run
    hid_t     file_id, dset_id, memspace, filespace, plist_id;
    hid_t     grp1,grp2,grp3,grp4;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel

#ifdef _MESSG
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
#endif
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //create file
    char fn_hd5[50];
    strcpy(fn_hd5,fbase);
    strcat(fn_hd5, ".h5");
    file_id = H5Fcreate(fn_hd5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);


    //1st 2 groups always active even if data is empty

    //group 1 is anything with mesh
    grp1 = H5Gcreate(file_id, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //group 2 contains any info with data
    grp2 = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //these groups on the other hand will always be there so this can be hard coded
    grp3 = H5Gcreate(grp1, "XY", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    grp4 = H5Gcreate(grp1, "Elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Gclose(grp1);
    H5Gclose(grp2);
    H5Gclose(grp3);
    H5Gclose(grp4);
    H5Fclose(file_id);
    
#endif

}

