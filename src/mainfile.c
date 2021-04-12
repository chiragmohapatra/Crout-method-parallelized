#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../include/matrix.h"

int n;  // square matrix of size n*n
double** A; // the input matrix
int num_threads;
int strategy;

void strat_0(double **A, double **L, double **U, int n) {
    int i, j, k;
    double sum = 0;
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            sum = 0;
            for (k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        for (i = j; i < n; i++) {
            sum = 0;
            for(k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0) {                
                exit(0);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
    }
}


int main(int argc , char* argv[]){
    if(argc != 5)
        exit(EXIT_FAILURE);

    n = atoi(argv[1]);

    FILE* fp = fopen(argv[2],"r");

    if(fp == NULL)
        exit(EXIT_FAILURE);

    char* line = NULL;
    size_t len = 0;
    A = allocate_matrix(n);

    int k = 0;
    
    while(getline(&line, &len, fp) != -1){
        assert(k < n);

        int l = 0;
        char* pch = strtok (line," ");

        while(pch != NULL){
            assert(l < n);

            A[k][l] = atof(pch);
            pch = strtok (NULL, " ");
            l++;
        }
        k++;
    }

    fclose(fp);
    if(line)
        free(line);

    num_threads = atoi(argv[3]);
    strategy = atoi(argv[4]);

    //print_matrix(A,n);

    double** L = allocate_matrix(n);
    double** U = allocate_matrix(n);

    if(strategy == 0){
        strat_0(A,L,U,n);
        print_matrix(L,n);
        print_matrix(U,n);
    }
    
    // write output to file
    char fname[] = {'o','u','t','p','u','t','_','L','_',(char)('0'+strategy),'_',(char)('0'+num_threads),'.','t','x','t','\0'};
    write_output(fname, L, n);
    fname[7] = 'U';
    write_output(fname, U, n);

    destroy_matrix(A,n);
    destroy_matrix(L,n);
    destroy_matrix(U,n);

    return 0;
}