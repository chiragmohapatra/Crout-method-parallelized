#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

// this functions allocates the memory for a square double matrix of size n*n
double** allocate_matrix(int n);

// this function deallocates the memory for a square matrix of size n*n
void destroy_matrix(double** A , int n);

// prints matrix in pretty format
void print_matrix(double** A , int n);

// write matrix contents to file
void write_output(char fname[], double** arr, int n );

#endif