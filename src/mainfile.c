#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

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

void strat_1(double **A, double **L, double **U, int n) {
    // setting up openmp
    omp_set_num_threads(num_threads);
    omp_set_nested(1); // for parallelising sum

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    // can't parallelize this outer loop since L[i][j] depends on L[i][k], k<j (and others)
    for (int j = 0; j < n; j++) {
        // static schedule since all iterations equally computationally expensive
        #pragma omp parallel for schedule(static)
        for (int i = j; i < n; i++) {
            double sum = 0;
            // parallelizing sum according to first strategy of assignment 1
            // #pragma omp parallel
            // {
            //     double partial_sum = 0;
            //     #pragma omp for
            //     for (int k = 0; k < j; k++) {
            //         partial_sum = partial_sum + L[i][k] * U[k][j];    
            //     }
            //     #pragma omp critical
            //     {
            //         sum += partial_sum;
            //     }
            // }
            for (int k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        #pragma omp parallel for schedule(static)
        for (int i = j; i < n; i++) {
            double sum = 0;
            // #pragma omp parallel
            // {
            //     double partial_sum = 0;
            //     #pragma omp for
            //     for (int k = 0; k < j; k++) {
            //         partial_sum = partial_sum + L[j][k] * U[k][i];    
            //     }
            //     #pragma omp critical
            //     {
            //         sum += partial_sum;
            //     }
            // }
            for(int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0) {                
                exit(0);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
    }
}

void strat_2(double **A, double **L, double **U, int n) {
    // setting up openmp
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        for (int i = 0; i < n; i++) {
            #pragma omp sections nowait
            {
                #pragma omp section
                U[i][i] = 1;
            }
        }
        #pragma omp barrier
    }

    for (int j = 0; j < n; j++) {
        #pragma omp parallel
        {
            for (int i = j; i < n; i++) {
                #pragma omp sections nowait
                {
                    #pragma omp section
                    {
                        double sum = 0;
                        for (int k = 0; k < j; k++) {
                            sum = sum + L[i][k] * U[k][j];    
                        }
                        L[i][j] = A[i][j] - sum;
                    }
                }
            }
            #pragma omp barrier

            for (int i = j; i < n; i++) {
                #pragma omp sections nowait
                {
                    #pragma omp section
                    {
                        double sum = 0;
                        for(int k = 0; k < j; k++) {
                            sum = sum + L[j][k] * U[k][i];
                        }
                        if (L[j][j] == 0) {                
                            exit(0);
                        }
                        U[j][i] = (A[j][i] - sum) / L[j][j];
                    }
                }
            }
            #pragma omp barrier
        }
    }
}

void strat_3(double **A, double **L, double **U, int n) {
    // setting up openmp
    omp_set_num_threads(num_threads);
    omp_set_nested(1);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    for (int j = 0; j < n; j++) {
        // calculate L[j][j] separately
        double sum = 0;
        #pragma omp parallel
        {
            double partial_sum = 0;
            #pragma omp for
            for (int k = 0; k < j; k++) {
                partial_sum = partial_sum + L[j][k] * U[k][j];    
            }
            #pragma omp critical
            {
                sum += partial_sum;
            }
        }
        L[j][j] = A[j][j] - sum;

        // calculate L[j..n][j] and U[j][j..n] in simultaneous sections 
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                #pragma omp parallel for schedule(static)
                for (int i = j; i < n; i++) {
                    double sum = 0;
                    // #pragma omp parallel
                    // {
                    //     double partial_sum = 0;
                    //     #pragma omp for
                    //     for (int k = 0; k < j; k++) {
                    //         partial_sum = partial_sum + L[i][k] * U[k][j];    
                    //     }
                    //     #pragma omp critical
                    //     {
                    //         sum += partial_sum;
                    //     }
                    // }
                    for (int k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];    
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }
            #pragma omp section
            {
                #pragma omp parallel for schedule(static)
                for (int i = j; i < n; i++) {
                    double sum = 0;
                    // #pragma omp parallel
                    // {
                    //     double partial_sum = 0;
                    //     #pragma omp for
                    //     for (int k = 0; k < j; k++) {
                    //         partial_sum = partial_sum + L[j][k] * U[k][i];    
                    //     }
                    //     #pragma omp critical
                    //     {
                    //         sum += partial_sum;
                    //     }
                    // }
                    for(int k = 0; k < j; k++) {
                        sum = sum + L[j][k] * U[k][i];
                    }
                    if (L[j][j] == 0) {                
                        exit(0);
                    }
                    U[j][i] = (A[j][i] - sum) / L[j][j];
                }
            }
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
    
    while(getline(&line, &len, fp) != -1 && k<n){
        assert(k < n);

        int l = 0;
        char* pch = strtok (line," ");

        while(pch != NULL && l<n){
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

    switch(strategy){
        case 0:
            strat_0(A,L,U,n);
            break;
        case 1:
            strat_1(A,L,U,n);
            break;
        case 2:
            strat_2(A,L,U,n);
            break;
        case 3:
            strat_3(A,L,U,n);
            break;
    }
    
    // print_matrix(L,n);
    // printf("\n");
    // print_matrix(U,n);

    // write output to file
    char fname[100] = {'o','u','t','p','u','t','_','L','_',(char)('0'+strategy),'_', '\0'};
    strcat(fname, argv[3]);
    strcat(fname, ".txt");
    write_output(fname, L, n);
    fname[7] = 'U';
    write_output(fname, U, n);

    destroy_matrix(A,n);
    destroy_matrix(L,n);
    destroy_matrix(U,n);

    return 0;
}