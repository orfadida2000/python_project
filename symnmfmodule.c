#include "symnmfmodule.h"

/* Calculates the symmetric matrix A and returns it as a Python list.
 * Uses `init_sym_matrix` to initialize A, and frees all allocated memory upon error.
 * Parameters:
 *  - self: Self-reference to the module (unused).
 *  - args: Python arguments passed from the caller.
 * Returns:
 *  - A_pylist: Python list representing the symmetric matrix A, or NULL on error.
 */
static PyObject* sym(PyObject* self, PyObject* args) {
    /* setup */
    PyObject *A_pylist;
    int n=0, flag = 0;
    double **n_mat=NULL, **A=NULL;
    /* inits the symmetrical matrix A */
    A = init_sym_matrix(args, &n_mat, &n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, NULL, NULL, n);
        return NULL;
    }
    /* transfrom A into a python list, returns NULL in case of error */
    A_pylist = mat_to_pylist(A, n, n);
    /* free all allocated memory using function from symnmf.c */
    free_all(n_mat, A, NULL, NULL, n);
    return A_pylist;

}

/* Calculates the diagonal degree matrix D and returns it as a Python list.
 * Initializes the symmetric matrix A and diagonal matrix D, and frees memory upon error.
 * Parameters:
 *  - self: Self-reference to the module (unused).
 *  - args: Python arguments passed from the caller.
 * Returns:
 *  - D_pylist: Python list representing the diagonal matrix D, or NULL on error.
 */
static PyObject* ddg(PyObject* self, PyObject* args) {
    /* setup */
    PyObject *D_pylist;
    int n=0, flag = 0;
    double **n_mat=NULL, **A=NULL, **D=NULL;
    /* inits the symmetrical matrix A */
    A = init_sym_matrix(args, &n_mat, &n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, D, NULL, n);
        return NULL;
    }
    D = calc_ddg_matrix(A, n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, D, NULL, n);
        return NULL;
    }
    /* transfrom  into a python list, returns NULL in case of error */
    D_pylist = mat_to_pylist(D, n, n);
    /* free all allocated memory using function from symnmf.c */
    free_all(n_mat, A, D, NULL, n);
    return D_pylist;

}

/* Calculates the normalized similarity matrix W and returns it as a Python list.
 * Initializes the symmetric matrix A, diagonal matrix D, and normalized matrix W.
 * Frees all allocated memory upon error.
 * Parameters:
 *  - self: Self-reference to the module (unused).
 *  - args: Python arguments passed from the caller.
 * Returns:
 *  - W_pylist: Python list representing the normalized matrix W, or NULL on error.
 */
static PyObject* norm(PyObject* self, PyObject* args) {
        /* setup */
    PyObject *W_pylist;
    int n=0, flag = 0;
    double **n_mat=NULL, **A=NULL, **D=NULL, **W=NULL;
    /* inits the symmetrical matrix A */
    A = init_sym_matrix(args, &n_mat, &n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, D, W, n);
        return NULL;
    }
    D = calc_ddg_matrix(A, n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, D, W, n);
        return NULL;
    }
    W = calc_norm_matrix(A, D, n, &flag);
    if (flag == 1) {
        free_all(n_mat, A, D, W, n);
        return NULL;
    }
    /* transfrom  into a python list, returns NULL in case of error */
    W_pylist = mat_to_pylist(W, n, n);
    /* free all allocated memory using function from symnmf.c */
    free_all(n_mat, A, D, W, n);
    return W_pylist;
}

/* Computes the final H matrix using SymNMF and returns it as a Python list.
 * Parses H and W from Python, allocates and computes final H matrix.
 * Parameters:
 *  - self: Self-reference to the module (unused).
 *  - args: Python arguments passed from the caller.
 * Returns:
 *  - H_pylist: Python list representing the computed H matrix, or NULL on error.
 */
static PyObject* symnmf(PyObject* self, PyObject* args) {
    PyObject *py_H =NULL, *py_W = NULL, *H_pylist = NULL;
    double **H0 = NULL, **H1 = NULL, **final_H = NULL, **W = NULL;
    double beta, eps; /* eps = epsilon */
    int n = 0, k = 0, flag = 0, maxIter;
    /* parse H and W from python into Py objects */
    py_H_W_parsing(args, &py_H, &py_W, &maxIter, &beta, &eps ,&flag);
    if (flag == 1)
        return NULL;
    /* transform H and W from PyObject* to double**, and assign n and k to correct values */
    H0 = py_Mat_To_Double_Mat(py_H, &n, &k, &flag); /* Initial H */
    W = py_Mat_To_Double_Mat(py_W, &n, &n, &flag);
    if (flag == 1) {
        free_all(H0, H1, W, NULL, n);
        return NULL;
    }
    /* setup for calculating final H */
    H1 = allocate_mat(n, k, &flag); /* allocates memory for n by k matrix (double**) and returns a pointer to it */
    if (flag == 1) {
        free_all(H0, H1, W, NULL, n);
        return NULL;
    }
    /* calculate final H and assigns it to either H0, or H1 */
    final_H = calc_H_mat(H0, H1, W, n, k, beta, eps, maxIter);
    /* transfrom into a python list, returns NULL in case of error */
    H_pylist = mat_to_pylist(final_H, n, k);
    /* free all allocated memory using function from symnmf.c */
    free_all(H0, H1, W, NULL, n);
    return H_pylist;
}

/* Calculates final H matrix using iterative updates until convergence or max iterations.
 * Parameters:
 *  - curr_H: Current H matrix (used in internal loop).
 *  - next_H: Next H matrix to compute (used in internal loop).
 *  - W: Similarity matrix.
 *  - n: Number of rows.
 *  - k: Number of columns.
 *  - beta: Update parameter.
 *  - eps: Convergence threshold.
 *  - max_iter: Maximum number of iterations.
 * Returns:
 *  - Final H matrix.
 */
static double** calc_H_mat(double** curr_H, double **next_H, double **W, int n, int k, double beta, double eps, int max_iter) {
    double **temp = NULL;
    int iter=1;
    Boolean conv = FALSE; /* We defined enum Boolean in symnmfmodule.h */
    while (conv == FALSE && iter <= max_iter){
        /* calculate next_H */
        calc_next_H(curr_H, next_H, W, beta, n, k);
        /* check convergence */
        conv = check_conv(curr_H, next_H, n, k, eps);
        /* update curr_H, and give next_H a different address */
        temp = curr_H;
        curr_H = next_H;
        next_H=temp;
        /* update iteration number */
        iter++;
    }
    return curr_H;
}

/* Checks convergence by comparing Frobenius norm difference of two matrices.
 * Parameters:
 *  - curr_H: Current H matrix.
 *  - next_H: Next H matrix.
 *  - rows, cols: Matrix dimensions.
 *  - eps: Convergence threshold.
 * Returns:
 *  - TRUE if converged, FALSE otherwise.
 */
static Boolean check_conv(double** curr_H, double **next_H, int rows, int cols, double eps) {
    if ( frobenius_metric_squared(curr_H, next_H, rows, cols) < eps ) {
        return TRUE;
    }
    return FALSE;
}

/* Calculates the squared Frobenius norm between two matrices A and B.
 * Parameters:
 *  - A, B: Matrices to compare.
 *  - rows, cols: Dimensions of matrices.
 * Returns:
 *  - Squared Frobenius norm as a double.
 */
static double frobenius_metric_squared(double **A, double **B, int rows, int cols) {
    int i,j;
    double res=0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            res += pow( A[i][j] - B[i][j] , 2 );
        }
    }
    return res;
}

/* Computes the next H matrix in the iterative SymNMF algorithm.
 * Parameters:
 *  - curr_H: Current H matrix.
 *  - next_H: Next H matrix to compute.
 *  - W: Normalized Similarity matrix.
 *  - beta: Update parameter.
 *  - rows, cols: Dimensions of H matrix.
 */
static void calc_next_H(double **curr_H, double **next_H, double **W, double beta, int rows, int cols) {
    int i,j,n = rows, k = cols;
    double numer, denom, entry_res;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            numer = calc_mat_mul_entry(W, curr_H, i, j, n, FALSE);
            denom = calc_H_denominator(curr_H,i,j,n,k);
            entry_res = curr_H[i][j] * (1-beta + beta*(numer/denom));
            next_H[i][j] = entry_res;
        }
    }
}

/* Calculates denominator for updating H matrix in SymNMF.
 * Parameters:
 *  - H: Current H matrix.
 *  - i, j: Indices in the matrix.
 *  - n, k: Matrix dimensions.
 * Returns:
 *  - Denominator value for updating H matrix.
 */
static double calc_H_denominator(double **H, int i, int j, int n, int k) {
    double res = 0, H_mul_Ht_entry_ip = 0; /* res, (H*H^t)[i,p] */
    int p;
    for (p=0;p<n;p++) {
        H_mul_Ht_entry_ip = calc_mat_mul_entry(H,H,i,p,k,TRUE); /* calculates (H*H^t)[i,p] */
        res += H_mul_Ht_entry_ip * H[p][j]; /* adds to the sum: ( (H*H^t)[i,p] * (H)[p,j] ) - in the sum over p */
    }
    return res;
}

/* Calculates the matrix multiplication entry of A and B (or A and B^t) for indices (i, j).
 * Parameters:
 *  - A, B: Matrices to multiply.
 *  - i, j: Indices for the entry.
 *  - mul_param: Dimension for multiplication.
 *  - is_B_Transposed: Boolean indicating if B is transposed - (calculates A * B^t if so).
 * Returns:
 *  - Result of matrix multiplication at entry (i, j).
 */
static double calc_mat_mul_entry(double **A, double **B, int i,int j, int mul_param, Boolean is_B_Transposed) {
    int A_row=i,A_col=0,B_row,B_col;
    double res=0;
    B_row = is_B_Transposed ? j : 0;
    B_col = is_B_Transposed ? 0 : j;
    while(A_col < mul_param) {
        res+= (A[A_row][A_col]*B[B_row][B_col]);
        if (is_B_Transposed){
            B_col++;
        } else{
            B_row++;
        }
        A_col++;
    }
    return res;
}

/* Allocates memory for a matrix of doubles with specified rows and columns.
 * Parameters:
 *  - rows, cols: Dimensions of the matrix.
 *  - flagPtr: Pointer to flag set to 1 on failure.
 * Returns:
 *  - Allocated matrix.
 */
static double** allocate_mat(int rows, int cols, int *flagPtr) {
    double** res = NULL;
    int i;
    res = (double**) calloc(rows, sizeof(double*));
    if (res == NULL) {
        *flagPtr = 1;
        return res;
    }
    for (i=0; i < rows; i++) {
        res[i] = malloc(cols * sizeof(double));
        if (res[i] == NULL) {
            *flagPtr = 1;
            return res;
        }
    }
    return res;
}

/* Verifies that a Python list is a list of lists.
 * Parameters:
 *  - py_List: Python object to check.
 *  - flagPtr: Pointer to flag set to 1 if not a list of lists, 0 otherwise.
 */
static void check_if_PyMatrix(PyObject* py_List, int* flagPtr) {
    PyObject *first_row = NULL;
    first_row = PyList_GetItem(py_List, 0);
    if (!PyList_Check(first_row)) { /* checks that the first row is a python list */
        *flagPtr = 1;
    }
}

/* Parses input Python objects and checks for errors.
 * Parameters:
 *  - args: Python arguments from the caller.
 *  - py_H_ptr, py_W_ptr: Pointers to Python lists representing H and W matrices.
 *  - iterPtr, betaPtr, epsPtr: Pointers to max iterations, beta, and epsilon values.
 *  - flagPtr: Pointer to flag set to 1 on failure.
 */
static void py_H_W_parsing(PyObject* args, PyObject** py_H_ptr,PyObject** py_W_ptr, int* iterPtr, double* betaPtr, double* epsPtr, int* flagPtr) {
    /* checks for errors in data received and assigns that data into the appropriate pointers */
    if (!PyArg_ParseTuple(args, "O!O!idd", &PyList_Type, py_H_ptr, &PyList_Type, py_W_ptr, iterPtr, betaPtr, epsPtr)) {
        *flagPtr = 1; /* raise error */
    }
    /* checks that H, and W are matrices (list of lists) as needed */
    check_if_PyMatrix(*py_H_ptr, flagPtr);
    check_if_PyMatrix(*py_W_ptr, flagPtr);
}

/* Parses input Python objects for the n_mat matrix and checks for errors.
 * Parameters:
 *  - args: Python arguments from the caller.
 *  - py_n_mat_ptr: Pointer to a Python list representing the n_mat matrix.
 *  - flagPtr: Pointer to a flag set to 1 on failure.
 */
static void py_nMat_parsing(PyObject* args, PyObject** py_n_mat_ptr, int* flagPtr) {
    /* checks for errors in data received and assigns that data into the appropriate pointers */
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, py_n_mat_ptr)) {
        *flagPtr = 1; /* raise error */
    }
    /* checks that n_mat (Data points matrix) is a matrix (list of lists) as needed */
    check_if_PyMatrix(*py_n_mat_ptr,flagPtr);
}

/* Initializes the symmetric matrix A from Python input.
 * Parameters:
 *  - args: Python arguments containing the n_mat matrix data.
 *  - n_mat_ptr: Pointer to the double matrix (n_mat).
 *  - nPtr: Pointer to the number of rows in n_mat.
 *  - flagPtr: Pointer to a flag set to 1 on error.
 * Returns:
 *  - A: Symmetric matrix as a double**.
 */
static double** init_sym_matrix(PyObject* args, double*** n_mat_ptr, int *nPtr, int *flagPtr) {
    PyObject *py_n_mat = NULL;
    double **A = NULL;
    int d;
    /* checks for errors in data received from python, and inits py_n_mat */
    py_nMat_parsing(args, &py_n_mat, flagPtr);
    if (*flagPtr == 1)
        return A;
    /* initializes basic data of data points */
    *n_mat_ptr = py_Mat_To_Double_Mat(py_n_mat, nPtr, &d, flagPtr);
    if (*flagPtr == 1)
        return A;

    /* calculates the symmetric matrix A */
    A = calc_sym_matrix(*n_mat_ptr, *nPtr, d, flagPtr);   /* calls function in symnmf.c, that calculates the symmetric matrix A */
    if (*flagPtr == 1)
        return A;

    return A;
}

/* Converts a Python matrix (list of lists) to a double matrix.
 * Parameters:
 *  - py_mat: Python matrix to convert.
 *  - rowsPtr, colsPtr: Pointers to store matrix dimensions.
 *  - flagPtr: Pointer to flag set to 1 on failure.
 * Returns:
 *  - Allocated double matrix.
 */
static double** py_Mat_To_Double_Mat(PyObject* py_mat, int *rowsPtr, int *colsPtr, int *flagPtr) {
    double **res = NULL;
    int i,j;
    *rowsPtr = (int) PyList_Size(py_mat);    /* number of data points */
    *colsPtr = (int) PyList_Size(PyList_GetItem(py_mat, 0));
    res = (double **) calloc(*rowsPtr, sizeof(double*));
    if (res == NULL) {
        *flagPtr = 1;
        return NULL;
    }

    for (i=0; i < *rowsPtr; i++) {
        res[i] = malloc(*colsPtr * sizeof(double));
        if (res[i] == NULL) {
            *flagPtr = 1;
            return res;
        }
        for (j=0; j < *colsPtr; j++) {
            res[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_mat, i), j));
        }
    }
    return res;
}

/* Converts a double matrix to a Python list of lists.
 * Parameters:
 *  - mat: Double matrix to convert.
 *  - rows, cols: Dimensions of the matrix.
 * Returns:
 *  - Python list of lists, or NULL on error.
 */
static PyObject* mat_to_pylist(double** mat, int rows, int cols) {
    if (mat == NULL) return NULL;
    int i = 0, j = 0;
    /* Create the outer Python list (new reference, refcount = 1) */
    PyObject* py_list = PyList_New(rows);
    if (!py_list) {
        return NULL;  /* Allocation failed, return an error */
    }
    /* Iterate over each row */
    for (i = 0; i < rows; i++) {
        /* Create a Python list for the current row (new reference, refcount = 1) */
        PyObject* py_row = PyList_New(cols);
        if (!py_row) {
            Py_DECREF(py_list);
            return NULL;  /* Allocation failed, return an error */
        }

        /* Populate the current row list with the elements from mat[i][] */
        for (j = 0; j < cols; j++) {
            PyObject* py_float = PyFloat_FromDouble(mat[i][j]);  /* New reference, refcount = 1 */
            if (!py_float) {
                Py_DECREF(py_row);
                Py_DECREF(py_list);
                return NULL;  /* Allocation failed, return an error */
            }
            PyList_SET_ITEM(py_row, j, py_float);  /* Steals reference to py_float, no need to Py_DECREF(py_float) */
        }
        /* Set the row in the outer list */
        PyList_SET_ITEM(py_list, i, py_row);  /* Steals reference to py_row, no need to Py_DECREF(py_row) */
    }
    /* Return the populated Python list of lists (refcount = 1, as expected) */
    return py_list;
}

/* Define the method table for the symnmf module */
static PyMethodDef SymnmfMethods[] = {

        {"sym",
         (PyCFunction) sym,
         METH_VARARGS,
         PyDoc_STR("Calculates the symmetric matrix A and returns it as a Python list.")},

         {"ddg",
         (PyCFunction) ddg,
         METH_VARARGS,
         PyDoc_STR("Calculates the diagonal degree matrix D and returns it as a Python list.")},

         {"norm",
         (PyCFunction) norm,
         METH_VARARGS,
         PyDoc_STR("Calculates the normalized similarity matrix W and returns it as a Python list.")},

         {"symnmf",
         (PyCFunction) symnmf,
         METH_VARARGS,
         PyDoc_STR("Computes the final H matrix using SymNMF and returns it as a Python list.")},

        {NULL, NULL, 0, NULL}  /* Sentinel */

};

/* Define the module */
static struct PyModuleDef symnmfmodule = {
        PyModuleDef_HEAD_INIT,
        "mysymnmfsp",
        NULL,
        -1,
        SymnmfMethods
};

/* Module initialization function */
PyMODINIT_FUNC PyInit_mysymnmfsp(void) {
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}
