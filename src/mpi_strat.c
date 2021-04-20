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

    double A[n][n];
    double L[n][n];
    double U[n][n];

    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            A[i][j] = 0;
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    // reading A from file
    if(my_rank == 0){

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
    }

    MPI_Bcast(&A , n*n , MPI_DOUBLE , 0 , MPI_COMM_WORLD);

    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    int elements_per_process = n/comm_sz;

    for (int j = 0; j < n; j++){
        double temp2[j]; // this is essentially U[k][j]
        double temp3[j+1]; // this is essentially L[j][k]

        if(my_rank == 0){
            // calculate L[j][j] early since this is used in U calculation
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][j];
            }
            L[j][j] = A[j][j] - sum;

            if(L[j][j] == 0){
                MPI_Finalize();
                exit(0);
            }

            for(int k = 0 ; k < j ; k++){
                temp2[k] = U[k][j];
                temp3[k] = L[j][k];
            }
            temp3[j] = L[j][j];
        }

        MPI_Bcast(&temp2, j, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&temp3, j+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(my_rank == 0){

            // now we will compute the part for process 0
            for(int i = j ; i < elements_per_process ; i++){
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

                int l = sender*elements_per_process;
                int u = (sender + 1)*elements_per_process;

                if(l < j)
                    l = j;

                if(u > n)
                    u = n;

                for (int i = l; i < u; i++){
                    if(tag == 0)
                        L[i][j] = sum_arr[i - l];
                    else
                        U[j][i] = sum_arr[i - l];
                }
            }
        }

        else{
            int l = my_rank*elements_per_process;
            int u = (my_rank + 1)*elements_per_process;

            if(l < j)
                l = j;

            if(u > n)
                u = n;

            int n_elements_received = u - l;
            if(u < l){
                double sum_arr1 = 0;
                MPI_Send(&sum_arr1,
                    1,MPI_DOUBLE,0,0,
                    MPI_COMM_WORLD);

                MPI_Send(&sum_arr1,
                    1,MPI_DOUBLE,0,1,
                    MPI_COMM_WORLD);
            }

            else{
                double sum_arr1[n_elements_received];

                for (int i = l; i < u; i++) {
                    double sum = 0;
                    for (int k = 0; k < j; k++) {
                        sum = sum + L[i][k] * temp2[k];    
                    }
                    L[i][j] = A[i][j] - sum;
                    sum_arr1[i - l] = L[i][j];
                }

                MPI_Send(&sum_arr1,
                    n_elements_received,MPI_DOUBLE,0,0,
                    MPI_COMM_WORLD);

                for (int i = l; i < u; i++) {
                    double sum = 0;
                    for (int k = 0; k < j; k++) {
                        sum = sum + U[k][i] * temp3[k];    
                    }
                    U[j][i] = (A[j][i] - sum) / temp3[j];
                    sum_arr1[i - l] = U[j][i];
                }

                MPI_Send(&sum_arr1,
                    n_elements_received,MPI_DOUBLE,0,1,
                    MPI_COMM_WORLD);
            }
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

        FILE* f;

        f = fopen(fname, "w");
        for( int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                fprintf(f, "%0.12f ", L[i][j]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        fname[7] = 'U';

        f = fopen(fname, "w");
        for( int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                fprintf(f, "%0.12f ", U[i][j]);
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    MPI_Finalize();
    return 0;
}