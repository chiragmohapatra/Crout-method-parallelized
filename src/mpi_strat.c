#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

// write matrix contents to file
void write_output(char fname[], double** arr, int n ){
	FILE *f = fopen(fname, "w");
	for( int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(f, "%0.12f ", arr[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

int main(int argc , char* argv[]){
    int n;  // square matrix of size n*n
    char* file_name; // name of the input file
    int num_threads;
    int strategy;

    if(argc != 4)
        exit(EXIT_FAILURE);

    n = atoi(argv[1]);
    file_name = argv[2];
    num_threads = atoi(argv[3]);
    strategy = 4;

    MPI_Status status;

    MPI_Init(NULL,NULL);

    int my_rank , comm_sz;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    double** A;
    double** L;
    double** U;

    // allocating memory for A,L,U and reading A from file
    if(my_rank == 0){
        A = malloc(n*sizeof(double*));
        L = malloc(n*sizeof(double*));
        U = malloc(n*sizeof(double*));

        for(int i = 0 ; i < n ; i++){
            A[i] = malloc(n*sizeof(double*));
            L[i] = malloc(n*sizeof(double*));
            U[i] = malloc(n*sizeof(double*));
        }

        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < n ; j++){
                A[i][j] = 0;
                L[i][j] = 0;
                U[i][j] = 0;
            }
        }

        FILE* fp = fopen(file_name,"r");

        if(fp == NULL)
            exit(EXIT_FAILURE);

        char* line = NULL;
        size_t len = 0;

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

        for (int i = 0; i < n; i++) {
            U[i][i] = 1;
        }
    }

    for (int j = 0; j < n; j++){
        double temp2[j]; // this is essentially U[k][j]
        double temp3[j]; // this is essentially L[j][k]

        if(my_rank == 0){
            // calculate L[j][j] early since this is used in U calculation
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][j];
            }
            L[j][j] = A[j][j] - sum;

            for(int k = 0 ; k < j ; k++){
                temp2[k] = U[k][j];
                temp3[k] = L[j][k];
            }
        }

        MPI_Bcast(&temp2, j, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&temp3, j, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(my_rank == 0){
            int elements_per_process = (n-j)/comm_sz;

            if(comm_sz > 1){
                for(int i = 1 ; i < comm_sz - 1 ; i++){
                    MPI_Send(&elements_per_process,
                         1, MPI_INT, i, 0,
                         MPI_COMM_WORLD);
                    
                    double temp1[elements_per_process][j]; // this is essentially L[i][k]
                    double temp4[j][elements_per_process]; // this is essentially U[k][i]
                    
                    int index = i*elements_per_process + j;

                    for(int k = index ; k < (index + elements_per_process) ; k++){
                        for(int l = 0 ; l < j ; l++){
                            temp1[k - index][l] = L[k][l];
                            temp4[l][k - index] = U[l][k];
                        }
                    }

                    MPI_Send(&temp1,
                        elements_per_process*j,MPI_DOUBLE,i,1,
                        MPI_COMM_WORLD);

                    MPI_Send(&temp4,
                        elements_per_process*j,MPI_DOUBLE,i,2,
                        MPI_COMM_WORLD);
                }

                // remaining for the last process
                int index = (comm_sz - 1)*elements_per_process + j;
                int elements_left = n - index;
                MPI_Send(&elements_left,
                    1, MPI_INT, comm_sz - 1, 0,
                    MPI_COMM_WORLD);

                double temp1[elements_left][j]; // this is essentially L[i][k]
                double temp4[j][elements_left]; // this is essentially U[k][i]

                for(int k = index ; k < n ; k++){
                    for(int l = 0 ; l < j ; l++){
                        temp1[k - index][l] = L[k][l];
                        temp4[l][k - index] = U[l][k];
                    }
                }

                MPI_Send(&temp1,
                    elements_left*j,MPI_DOUBLE,comm_sz - 1,1,
                    MPI_COMM_WORLD);

                MPI_Send(&temp4,
                    elements_left*j,MPI_DOUBLE,comm_sz - 1,2,
                    MPI_COMM_WORLD);
            }

            // now we will compute the part for process 0
            for(int i = j ; i < j + elements_per_process ; i++){
                double sum1 = 0 , sum2 = 0;
                for (int k = 0; k < j; k++) {
                    sum1 = sum1 + L[i][k] * U[k][j];  
                    sum2 = sum2 + L[j][k] * U[k][i];  
                }
                L[i][j] = A[i][j] - sum1;
                U[j][i] = (A[j][i] - sum2) / L[j][j];
            }

            // Now we will collect the computations from all the other processes
            for(int i = 0 ; i < 2*(comm_sz - 1) ; i++){
                double sum_arr[n];

                MPI_Recv(&sum_arr,
                    n,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
                    MPI_COMM_WORLD,
                    &status);

                int sender = status.MPI_SOURCE;
                int tag = status.MPI_TAG;

                int index = sender*elements_per_process + j;

                int u = index + elements_per_process;

                if(sender == comm_sz - 1)
                    u = n;

                for(int k = index ; k < u ; k++){
                    if(tag == 0)
                        L[k][j] = A[k][j] - sum_arr[k - index];
                    else{
                        U[j][k] = (A[j][k] - sum_arr[k - index]) / L[j][j];
                    }
                }
            }
        }

        else{
            int n_elements_recieved;

            MPI_Recv(&n_elements_recieved,
                 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD,
                 &status);

            double temp1[n_elements_recieved][j]; // this is essentially L[i][k]
            double temp4[j][n_elements_recieved]; // this is essentially U[k][i]

            MPI_Recv(&temp1,
                n_elements_recieved*j,MPI_DOUBLE,0,1,
                MPI_COMM_WORLD,
                &status);

            double sum_arr1[n_elements_recieved];

            for(int i = 0 ; i < n_elements_recieved ; i++){
                double sum1 = 0;
                for (int k = 0; k < j; k++) {
                    sum1 = sum1 + temp1[i][k] * temp2[k];    
                }
                sum_arr1[i] = sum1;
            }

            MPI_Send(&sum_arr1,
                n_elements_recieved,MPI_DOUBLE,0,0,
                MPI_COMM_WORLD);

            MPI_Recv(&temp4,
                n_elements_recieved*j,MPI_DOUBLE,0,2,
                MPI_COMM_WORLD,
                &status);

            //double sum_arr2[n_elements_recieved];

            for(int i = 0 ; i < n_elements_recieved ; i++){
                double sum1 = 0;
                for (int k = 0; k < j; k++) {
                    sum1 = sum1 + temp4[k][i] * temp3[k];    
                }
                sum_arr1[i] = sum1;
            }

            MPI_Send(&sum_arr1,
                n_elements_recieved,MPI_DOUBLE,0,1,
                MPI_COMM_WORLD);
        }
    }

    if(my_rank == 0){
        // write output to file
        char fname[100] = {'o','u','t','p','u','t','_','L','_',(char)('0'+strategy),'_'};
        char str_threads[5];
        sprintf(str_threads , "%d", num_threads);
        char temp[] = {'.','t','x','t','\0'};

        strcat(fname,str_threads);
        strcat(fname,temp);

        write_output(fname, L, n);
        fname[7] = 'U';
        write_output(fname, U, n);

        for(int i = 0 ; i < n ; i++){
            free(A[i]);
            free(L[i]);
            free(U[i]);
        }

        free(A);
        free(L);
        free(U);
    }

    MPI_Finalize();
    return 0;
}