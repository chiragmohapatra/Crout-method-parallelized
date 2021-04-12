#include "../include/matrix.h"

// this functions allocates the memory for a square double matrix of size n*n
double** allocate_matrix(int n){
    double** A = (double**)malloc(n*sizeof(double*));

    for(int i = 0 ; i < n ; i++){
        A[i] = (double*)malloc(n*sizeof(double));
    }

    return A;
}

// this function deallocates the memory for a square matrix of size n*n
void destroy_matrix(double** A , int n){
    for(int i = 0 ; i < n ; i++){
        free(A[i]);
    }

    free(A);
}

// prints matrix in pretty format
void print_matrix(double** A , int n){
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            printf("%lf " , A[i][j]);
        }
        printf("\n");
    }
}