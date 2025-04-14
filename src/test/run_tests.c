/*! \file  run_tests.c This file has calls all unit tests*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function calls all engine tests
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \returns 0 if succesful, nonzero if something failed
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int run_tests(void){

    int all_err = 0;
    int err = 0;

if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Beggining Unit Tests\n");
        printf("------------------------------------------------------\n");
}

#ifdef _MESSG
        all_err = test_comm_update(P2P);
        all_err += test_comm_update(NEIGHBOR);
        all_err += test_grid_partition(NEIGHBOR);
#else
if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("Residual Tests Begin\n");
    printf("------------------------------------------------------\n");
    
}
    int n_resid_tests = 1;
    int npx[] = {3,20,4,5,10,11,8,6,9,350};
    int npy[] = {3,20,4,3,12,8,9,7,6,300};
    double xmin[] = {0.0, -2.0, -5.0, 1.0, -2.0, -10.0, 2.0, -100.0, 29.0, 21.0 };
    double xmax[] = {1.0, 1.0, 2.0, 3.0, 0.0, -5.0, 4.0, -98.5, 32.0, 24.0};
    double ymin[] = {0.0, -2.0, -5.0, 1.2, -10.0, 25.0, 100.0, -1000.0, 0.0, 1.25};
    double ymax[] = {1.0, 1.0, 2.0, 1.8, -8.5, 26.3, 102.5, -998.5, 5.0, 3.75};
    int nts[] = {50, 20, 10, 15, 16, 30, 7, 18, 13, 20};
    for (int i = 0 ; i < n_resid_tests; i++){
        err += test_residual(npx[i],npy[i],xmin[i],xmax[i],ymin[i],ymax[i]);
    }
   if(err!=0){all_err+=1;}
if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("%d / %d Residual Tests Passed\n", err+n_resid_tests, n_resid_tests);
    printf("------------------------------------------------------\n"); 
}
if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("Jacobian Tests Begin\n");
    printf("------------------------------------------------------\n"); 
}
    err = 0;
    int n_jacobian_tests = 10;//
    for (int i = 0 ; i < n_jacobian_tests; i++){
        err += test_jacobian(npx[i],npy[i],xmin[i],xmax[i],ymin[i],ymax[i]);//
    }
    if(err!=0){all_err+=1;}
if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("%d / %d Jacobian Tests Passed\n", err+n_jacobian_tests, n_jacobian_tests);
    printf("------------------------------------------------------\n");   
}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("Newton Tests Begin\n");
//    printf("------------------------------------------------------\n");  
//}
//    err = 0;
//    int n_newton_tests = 10;
//    for (int i = 0 ; i < n_newton_tests; i++){
//        err += test_newton(npx[i],npy[i],xmin[0],xmax[0],ymin[0],ymax[0]);
//        printf(" NEWTON TEST %d / %d completed \n",i+1,n_newton_tests);
//    }
//    if(err!=0){all_err+=1;}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("%d / %d Newton Tests Passed\n", err+n_jacobian_tests, n_jacobian_tests);
//    printf("------------------------------------------------------\n");
//    
//}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("Time Loop Tests Begin\n");
//    printf("------------------------------------------------------\n"); 
//}
//    err = 0;
//    int n_timeloop_tests = 10;
//    npx[9] = 300;
//    npy[9] = 350;  
//    for (int i = 0 ; i < n_timeloop_tests; i++){
//        err += test_timeloop(npx[i],npy[i], nts[i]);
//        printf(" Timeloop TEST %d / %d completed \n",i+1,n_timeloop_tests);
//    }
//    if(err!=0){all_err+=1;}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("%d / %d Timeloop Tests Passed\n", err+n_timeloop_tests, n_timeloop_tests);
//    printf("------------------------------------------------------\n");   
//}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("SW2 WD Tests Begin\n");
//    printf("------------------------------------------------------\n"); 
//}
//    err = 0;
//    int n_sw2_wd_tests = 1;
//    npx[0] = 16;
//    npy[0] = 6;
//    nts[0] = 1440;
//    for (int i = 0 ; i < n_sw2_wd_tests; i++){
//        err += test_sw2_wd(npx[i],npy[i], nts[i]);
//        printf(" sw2_wd_test TEST %d / %d completed \n",i+1,n_sw2_wd_tests);
//    }
//    if(err!=0){all_err+=1;}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("%d / %d Sw2 wd Tests Passed\n", err+n_sw2_wd_tests, n_sw2_wd_tests);
//    printf("------------------------------------------------------\n");
//    
//}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("SW2 NB Tests Begin\n");
//    printf("------------------------------------------------------\n");    
//}
//    err = 0;
//    int n_sw2_nb_tests = 1;
//    nts[0] = 10;
//    for (int i = 0 ; i < n_sw2_nb_tests; i++){
//        err += test_sw2_nb(nts[i]);
//        printf(" sw2_nb_test TEST %d / %d completed \n",i+1,n_sw2_nb_tests);
//    }
//    if(err!=0){all_err+=1;}
//if (DEBUG) {
//    printf("------------------------------------------------------\n");
//    printf("------------------------------------------------------\n");
//    printf("%d / %d Sw2 nb Tests Passed\n", err+n_sw2_nb_tests, n_sw2_nb_tests);
//    printf("------------------------------------------------------\n");
//    
//}
#endif





if (DEBUG) {
    if (all_err == 0){
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("All Unit Tests Passed\n");
        printf("------------------------------------------------------\n");
    }else{
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Some Unit Tests Failed\n");
        printf("------------------------------------------------------\n");
    }
    
}


	return all_err; 

}