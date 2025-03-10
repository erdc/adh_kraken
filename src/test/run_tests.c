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

    int err = 0;

if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Beggining Unit Tests\n");
        printf("------------------------------------------------------\n");
}

if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("Residual Tests Begin\n");
    printf("------------------------------------------------------\n");
    
}
    int n_resid_tests = 3;
    int npx[] = {3,20,4};
    double xmin[] = {0.0, -2.0, -5.0};
    double xmax[] = {1.0, 1.0, 2.0};

    for (int i = 0 ; i < n_resid_tests; i++){
        err += residual_test(npx[i],xmin[i],xmax[i]);
    }


if (DEBUG) {
    printf("------------------------------------------------------\n");
    printf("------------------------------------------------------\n");
    printf("%d / %d Residual Tests Passed\n", err+n_resid_tests, n_resid_tests);
    printf("------------------------------------------------------\n");
    
}




	return err; 

}