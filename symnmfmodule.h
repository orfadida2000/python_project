#ifndef SYMNMF_MODULE_H
#define SYMNMF_MODULE_H

#include <Python.h>  /* Include the Python C API */
#include "symnmf.h"  /* Include the symnmf header file */

/* Define a Boolean enum for ANSI C */
typedef enum {
    FALSE = 0,
    TRUE = 1
} Boolean;

/* Function prototypes for functions exposed to Python */
static PyObject* sym(PyObject* self, PyObject* args);
static PyObject* ddg(PyObject* self, PyObject* args);
static PyObject* norm(PyObject* self, PyObject* args);
static PyObject* symnmf(PyObject* self, PyObject* args);

/* Data parsing from python */
static void py_nMat_parsing(PyObject* args, PyObject** py_n_mat_ptr, int* flagPtr);
static void py_H_W_parsing(PyObject* args, PyObject** py_H_ptr,PyObject** py_W_ptr, int* iterPtr, double* betaPtr, double* epsPtr, int* flagPtr);
static void check_if_PyMatrix(PyObject* py_List, int* flagPtr);
static double** py_Mat_To_Double_Mat(PyObject* py_mat, int *rowsPtr, int *colsPtr, int *flagPtr);
static PyObject* mat_to_pylist(double** mat, int rows, int cols);

/* symnmf functions */
static double** calc_H_mat(double** curr_H, double **next_H, double **W, int n, int k, double beta, double eps, int max_iter);
static double calc_H_denominator(double **H, int i, int j, int n, int k);
static void calc_next_H(double **curr_H, double **next_H, double **W, double beta, int rows, int cols);
static Boolean check_conv(double** curr_H, double **next_H, int rows, int cols, double eps);

/* General matrix functions */
static double** allocate_mat(int rows, int cols, int *flagPtr);
static double calc_mat_mul_entry(double **A, double **B, int i,int j, int mul_param, Boolean is_B_Transposed);
static double frobenius_metric_squared(double **A, double **B, int rows, int cols);

/* the symmetrical matrix A */
static double** init_sym_matrix(PyObject* args, double*** n_mat_ptr, int *nPtr, int *flagPtr);

#endif  /* SYMNMF_MODULE_H */
