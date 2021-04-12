#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../include/matrix.h"

int n;  // square matrix of size n*n
double** A; // the input matrix
int num_threads;
int strategy;


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

    print_matrix(A,n);

    destroy_matrix(A,n);

    return 0;
}