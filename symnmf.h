#ifndef SYM_NMF_H
#define SYM_NMF_H

#ifdef _WIN32
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>  /* Standard library for memory allocation */
#include <stdio.h>   /* Standard I/O */
#include <math.h>    /* Math library */
#include <string.h>  /* String library */

/* init n_mat helper functions */
int gline(char **lineptr, size_t *n, FILE *fp);
int init_nmat(double ***n_matPtr, int *nptr, int *dptr, FILE *fp);
int get_d(const char *line);
void line_to_vector(char *line, double *vector, int d);
int get_file_data(const char *file_name,double ***n_matPtr,int *nptr, int *dptr);
void handle_error_and_exit(double **n_mat, double **A, double **D, double **W, int n);

/* output or ending functions */
void print_matrix(double **matrix, int rows, int cols);
void free_matrix(double **matrix, int rows);
void free_all(double **n_mat, double **A, double **D, double **W, int n);

/* general helper function */
double distance_squared(double *u, double *v, int d);

/* symNMF calculating functions */
double** calc_sym_matrix(double **n_mat, int n, int d, int *flag);
double** calc_ddg_matrix(double **A, int n, int *flag);
double** calc_norm_matrix(double** A, double** D,int n, int *flag);

#endif