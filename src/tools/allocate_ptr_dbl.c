#include "adh.h"

int **allocate_dptr_int(int rows, int cols) {

    char *str = NULL; // for error messages
    int **arr = (int**) tl_alloc(sizeof(int *), rows); // Allocate memory for row pointers
    if (arr == NULL) {
        sprintf(str, "Memory Allocation Failure\n");
        tl_error(str);
    }

    for (int i = 0; i < rows; i++) {
        arr[i] = (int*) tl_alloc(sizeof(int), cols); // Allocate memory for each row

        if (arr[i] == NULL) {
            // If memory allocation fails for a row, free previously allocated memory
            for (int j = 0; j < i; j++) {
                arr[j] = (int *) tl_free(sizeof(int),cols,arr[j]);
            }
            arr = (int **) tl_free(sizeof(int *),rows,arr);
            sprintf(str, "Memory Allocation Failure\n");
            tl_error(str);
        }
    }

    return arr;
}
void free_dptr_int(int **ptr, int rows, int cols) {
    if (ptr == NULL) return; // Check for NULL pointer

    // Free each row
    for (int i = 0; i < rows; i++) {
        ptr[i] = (int *) tl_free(sizeof(int),cols,ptr[i]);
    }

    // Free the array of pointers
    ptr = (int **) tl_free(sizeof(int *),rows,ptr);
}

char **allocate_dptr_char(int rows, int cols) {

    char *str = NULL; // for error messages
    char **arr = (char**) tl_alloc(sizeof(char *), rows); // Allocate memory for row pointers
    if (arr == NULL) {
        sprintf(str, "Memory Allocation Failure\n");
        tl_error(str);
    }

    for (int i = 0; i < rows; i++) {
        arr[i] = (char*) tl_alloc(sizeof(char), cols); // Allocate memory for each row

        if (arr[i] == NULL) {
            // If memory allocation fails for a row, free previously allocated memory
            for (int j = 0; j < i; j++) {
                arr[j] = (char *) tl_free(sizeof(char),cols,arr[j]);
            }
            arr = (char **) tl_free(sizeof(char *),rows,arr);
            sprintf(str, "Memory Allocation Failure\n");
            tl_error(str);
        }
    }

    return arr;
}
void free_dptr_char(char **ptr, int rows, int cols) {
    if (ptr == NULL) return; // Check for NULL pointer

    // Free each row
    for (int i = 0; i < rows; i++) {
        ptr[i] = (char *) tl_free(sizeof(char),cols,ptr[i]);
    }

    // Free the array of pointers
    ptr = (char **) tl_free(sizeof(char *),rows,ptr);
}



double **allocate_dptr_dbl(int rows, int cols) {

    char *str = NULL; // for error messages
    double **arr = (double **) tl_alloc(sizeof(double *), rows); // Allocate memory for row pointers
    if (arr == NULL) {
        sprintf(str, "Memory Allocation Failure\n");
        tl_error(str);
    }

    for (int i = 0; i < rows; i++) {
        arr[i] = (double *) tl_alloc(sizeof(double), cols); // Allocate memory for each row

        if (arr[i] == NULL) {
            // If memory allocation fails for a row, free previously allocated memory
            for (int j = 0; j < i; j++) {
                arr[j] = (double *) tl_free(sizeof(double),cols,arr[j]);
            }
            arr = (double **) tl_free(sizeof(double *),rows,arr);
            sprintf(str, "Memory Allocation Failure\n");
            tl_error(str);
        }
    }

    return arr;
}
void free_dptr_dbl(double **ptr, int rows, int cols) {
    if (ptr == NULL) return; // Check for NULL pointer

    // Free each row
    for (int i = 0; i < rows; i++) {
        ptr[i] = (double *) tl_free(sizeof(double),cols,ptr[i]);
    }

    // Free the array of pointers
    ptr = (double **) tl_free(sizeof(double *),rows,ptr);
}
