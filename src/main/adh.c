/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This file containers the stand-alone driver to AdH
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The stand-alone driver for AdH
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] argc           (int)  # of command line arguements
 * @param[in] argv           (char **)  command line arguements
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static time_t time1, time2, time_start;        /* for run-time calculation */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/\
int main(int argc, char *argv[]) {
    
    bool input_check = false;
    int i, myid = 0, npes = 1, ierr_code = UNSET_INT;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize MPI
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _MESSG
    MPI_Comm comm_world = MPI_COMM_WORLD;
    ierr_code = MPI_Init(&argc, &argv);
    if (ierr_code != MPI_SUCCESS) {messg_err(ierr_code);}
    //ierr_code = MPI_Comm_rank(cstorm_comm, &myid); if (ierr_code != MPI_SUCCESS) {messg_err(ierr_code);}
    //ierr_code = MPI_Comm_size(cstorm_comm, &npes); if (ierr_code != MPI_SUCCESS) {messg_err(ierr_code);}
    ierr_code = MPI_Comm_rank(comm_world, &myid); if (ierr_code != MPI_SUCCESS) {messg_err(ierr_code);}
    ierr_code = MPI_Comm_size(comm_world, &npes); if (ierr_code != MPI_SUCCESS) {messg_err(ierr_code);}
#endif

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize the AdH memory debugger
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _MESSG
    debug_initialize(comm_world);
#else
    debug_initialize();
#endif
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize PETSC
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef _PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);
#endif
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize all function pointers
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    set_function_pointers(adh_resid_routines, adh_init_routines, adh_time_stepper);

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Check command line arguments
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        for (i=0; i<argc; i++) {printf("command line argument[%d]: %s\n",i,argv[i]);}
    }
    if (argc < 2 || argc > 5) {
        fprintf(stderr, "\n Missing argument.\n Usage:\n" "   To run a simulation:\n     adh file_base -s (optional)\n" "   Or,\n   For version information:\n      adh -v\n");
        exit(0);
        
    } else if (strcmp(argv[1], "-t") == AGREE) {
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++ TESTING ADH BUILD +++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        run_tests();
#ifdef _PETSC
    PetscFinalize();
#endif
#ifdef _MESSG
    MPI_Finalize();
#endif
        exit(0);
        
    } else if (strcmp(argv[1], "-v") == AGREE) {
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("++++++++++++++++ ADH VERSION AND BUILD DATA +++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        print_build_info();
        printf("\n");
        exit(0);
    } else if (strcmp(argv[1], "-s") == AGREE) {
        printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("++++++++++++++ CHECKING ADH INPUT FILES ONLY +++++++++++++++\n");
        printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        input_check = true;
    }
    
    char filename[30];
    strcpy(filename,argv[1]); printf("filename: %s\n",filename);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // allocate an AdH design model
    SMODEL_DESIGN dmod; smodel_design_defaults(&dmod);
    strcpy(dmod.filebase,filename);

    // start calculation time
    time(&time_start);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // AdH initialization
    time(&time1);
    int ierr = 0;
#ifdef _MESSG
    ierr = smodel_design_init(&dmod, filename, input_check, &comm_world);
#else
    ierr = smodel_design_init(&dmod, filename, input_check);
#endif
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize HDF5/XDMF file output
#ifdef _HDF5
    char domain_name[50] = "adh sim";
    smodel_design_xmf_init(&dmod, filename, domain_name);

    // Write initial time-step ata
    smodel_design_xmf_write(&dmod,0); // write intial data
    
    // TESTING -- write a second time-step
    //dmod.ts.nt = 1;
    //dmod.ts.time = 100.0;
    //smodel_design_xmf_write(&dmod,0);
#endif

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // End AdH initialization
    time(&time2);
    double time_initial = difftime(time2, time1);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    
//    // AdH execution
//    time(&time1);
//    ierr = adh_run_func_(superModel, superinterface, nsupermodels, nsuperinterfaces, &run_time);
//    time(&time2);
//    double time_run = difftime(time2, time1);
//
//    int ierr_code, myid=0;
//#ifdef _MESSG
//    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//#endif
//    if (myid <= 0) {
//        printf("\n");
//#ifdef _DEBUG
//        printf("**********************************************************************\n");
//        printf("Timings: \n");
//        printf("-- Total initialization time: %lf seconds\n", time_initial);
//        for (i=0; i<nsupermodels; i++) {
//            for (j=0; j<superModel[i].nsubmodels; j++) {
//                if (superModel[i].submodel[j].flag.SW2_FLOW) {
//                    printf("SW 2D TIMINGS\n");
//                    printf("---- HYD RESID runtime: %lf seconds\n", TIME_IN_2D_SW_RESID);
//                    printf("---- HYD LOAD runtime: %lf seconds\n", TIME_IN_2D_SW_LOAD);
//                    printf("------ HYD 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_2D_SW_BODY_RESID);
//                    printf("------ HYD 1D ELEM RESID runtime: %lf seconds\n", TIME_IN_2D_SW_BOUNDARY_RESID);
//                } else if (superModel[i].submodel[j].flag.SW3_FLOW) {
//                    printf("SW 3D TIMINGS\n");
//                    printf("---- HVEL RESID runtime: %lf seconds\n", TIME_IN_HVEL_RESID);
//                    printf("---- HVEL LOAD runtime: %lf seconds\n", TIME_IN_HVEL_LOAD);
//                    printf("------ HVEL 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_HVEL_BODY_RESID);
//                    printf("------ HVEL 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_HVEL_BOUNDARY_RESID);
//                    printf("---- WVEL RESID runtime: %lf seconds\n", TIME_IN_WVEL_RESID);
//                    printf("---- WVEL LOAD runtime: %lf seconds\n", TIME_IN_WVEL_LOAD);
//                    printf("------ WVEL 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_WVEL_BODY_RESID);
//                    printf("------ WVEL 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_WVEL_BOUNDARY_RESID);
//                } else if (superModel[i].submodel[j].flag.GW_FLOW) {
//                    printf("---- GW RESID runtime: %lf seconds\n", TIME_IN_GW_RESID);
//                    printf("---- GW LOAD runtime: %lf seconds\n", TIME_IN_GW_LOAD);
//                    printf("------ GW 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_GW_BODY_RESID);
//                    printf("------ GW 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_GW_BOUNDARY_RESID);
//
//                }
//            }
//        }
//#endif
//    }
//
//    // AdH finalization
//    time(&time1);
//    ierr = adh_finalize_func_(&superModel, &superinterface, nsupermodels, nsuperinterfaces);
//    time(&time2);
//    double time_final = difftime(time2, time1);
//
//    /*******************************************************/
//    /*******************************************************/
//    /* getting calculation time */
//    time(&time2);
//    double t2 = difftime(time2, time_start);
//
//    if (myid <= 0) {
//#ifdef _DEBUG
//        printf("-- Total execution time: %lf seconds\n", time_run);
//        printf("-- Total finalization time: %lf seconds\n", time_final);
//#endif
//        printf("Total simulation runtime is %lf seconds\n", t2);
//    }

//#ifdef _PETSC
//    ierr = PetscFinalize();
//#endif
//#ifdef _MESSG
//    MPI_Finalize();
//#endif
//
#ifdef _HDF5
    xmf_finalize(dmod.xmf, dmod.grid->smpi->myid);
#endif
#ifdef _PETSC
    PetscFinalize();
#endif

    return ierr;
}

/*******************************************************/
/*******************************************************/
/*******************************************************/
