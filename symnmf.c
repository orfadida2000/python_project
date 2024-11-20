#include "symnmf.h"

/* Main function to execute the program.
 * Takes in command-line arguments specifying the goal and file name.
 * Valid goals: "sym", "ddg", "norm"
 * Returns 0 on success, exits with error code 1 on failure.
 */
int main(int argc, char *argv[]) {
    int n=0, d, flag = 0;  /* n = number of vectors, d = dimension of vectors */
    double **n_mat = NULL, **A = NULL, **D = NULL, **W = NULL;
    char *goal = argv[1];
    if(argc!=3){
        handle_error_and_exit(n_mat, A, D, W, n);
    }
    if ((get_file_data(argv[2], &n_mat, &n, &d)) == -1) {
        handle_error_and_exit(n_mat, A, D, W, n);
    }
    if ((A=calc_sym_matrix(n_mat, n, d,&flag)),flag == 1) {
        handle_error_and_exit(n_mat, A, D, W, n);
    }
    if (strcmp(goal, "sym") == 0) {
        print_matrix(A, n, n);
    }
    else if (strcmp(goal, "ddg") == 0) {
        if ((D=calc_ddg_matrix(A, n,&flag)),flag == 1) {
            handle_error_and_exit(n_mat, A, D, W, n);
        }
        print_matrix(D, n, n);
    }
    else if (strcmp(goal, "norm") == 0) {
        if ( (D=calc_ddg_matrix(A, n,&flag)),flag == 1 || (W = calc_norm_matrix(A, D, n,&flag)),flag == 1) {
            handle_error_and_exit(n_mat, A, D, W, n);
        }
        print_matrix(W, n, n);
    }
    free_all(n_mat, A, D, W, n);
    return 0;
}

/* Reads data from a file to populate the n_mat matrix with vectors.
 * Arguments:
 * file_name : name of the file to read from
 * n_matPtr  : pointer to the matrix that stores the vectors
 * nptr      : pointer to number of stored vectors
 * dptr      : pointer to stored vectors' dimension
 * Returns 1 on success, -1 on failure.
 */
int get_file_data(const char *file_name, double ***n_matPtr, int *nptr, int *dptr) {
    FILE *file = fopen(file_name, "r");
    if (file == NULL)
        return -1;
    if (init_nmat(n_matPtr, nptr, dptr, file) == -1)
        return -1;
    fclose(file);
    return 1;
}

/* Initializes all input vectors into a dynamically allocated matrix n_mat.
 * Returns 0 on success, -1 on failure.*/
int init_nmat(double ***n_matPtr, int *nptr, int *dptr, FILE *fp) {
    double **final_n_mat = NULL;
    int count = 0, read, d = -1, n = 5;
    char *line = NULL;
    size_t len = 0;
    if (n_matPtr == NULL) return -1;
    if (*n_matPtr == NULL) { /* Initial allocation */
        *n_matPtr = (double **) malloc(n * sizeof(double *));
        if (*n_matPtr == NULL) return -1;
    }
    while ((read = gline(&line, &len, fp)) != -1 && read != -2) /* gline can return negative number meaning error */
    {
        if (d == -1) d = get_d(line); /* Find d*/
        if (count >= n) { /* Doubling n_mat size */
            int new_size = n * 2;
            double **new_n_mat = (double **) realloc(*n_matPtr, new_size * sizeof(double *));
            if (new_n_mat == NULL) {
                *nptr = count;
                return -1;
            }
            *n_matPtr = new_n_mat;
            n = new_size;
        }
        (*n_matPtr)[count] = (double *) malloc(d * sizeof(double)); /* allocate vector memory */
        if ((*n_matPtr)[count] == NULL) {
            *nptr = count;
            return -1; /*Need to free everything before*/
        }
        line_to_vector(line, (*n_matPtr)[count], d);
        count++;
    }
    free(line);
    if (read == -1) return -1; /* error in gline */
    *nptr = count, *dptr = d; /* update global variables n , d */
    final_n_mat = (double **) realloc(*n_matPtr, count * sizeof(double *)); /* realloc n_mat to final size */
    if (final_n_mat == NULL) return -1;
    *n_matPtr = final_n_mat;
    return 0;
}

/* Determines the dimension 'd' by counting commas in the first line of the file.
 * Arguments:
 * line : a single line from the file
 * Returns the dimension d as an integer.
 */
int get_d(const char *line) {
    int i = 0, d = 1;
    while (line[i] != '\0') {
        if (line[i++] == ',')
            d++;
    }
    return d;
}

/* Converts a comma-separated line of numbers into a vector.
 * Arguments:
 * line   : string containing comma-separated values
 * vector : array to store parsed values
 * d      : number of elements in the vector
 */
void line_to_vector(char *line, double *vector, int d) {
    char *start = line;
    char *end = NULL;
    double curr_num;
    int i;
    for (i = 0; i < d; i++) {
        curr_num = strtod(start, &end);
        start = end + 1;
        vector[i] = curr_num;
    }
}

/* Reads a line from a file, dynamically resizing the buffer as needed.
 * Arguments:
 * lineptr : pointer to the buffer to store the line
 * n       : pointer to the size of the buffer
 * fp      : file pointer
 * Returns the number of characters read, or -1 on failure, or -2 when line is empty
 */
int gline(char **lineptr, size_t *n, FILE *fp) {
    size_t pos = 0, final_size;
    char *final_lineptr = NULL;
    int c;
    if (lineptr == NULL || n == NULL) return -1;
    if (*lineptr == NULL) {
        *n = 128;
        *lineptr = malloc(*n);
        if (*lineptr == NULL)
            return -1;
    }
    while ((c = fgetc(fp)) != '\n' && c != EOF) {
        /* double lineptr size */
        if (pos + 1 >= *n) {
            size_t new_size = *n * 2;
            char *new_lineptr = (char *) realloc(*lineptr, new_size);
            if (new_lineptr == NULL)
                return -1;
            *lineptr = new_lineptr;
            *n = new_size;
        }
        (*lineptr)[pos++] = (char) c;
    }
    if (pos == 0) /* check if enter here ever */
        return -2;
    (*lineptr)[pos] = '\0';
    final_size = pos + 1; /* Include space for the null terminator*/
    final_lineptr = (char *) realloc(*lineptr, final_size);
    if (final_lineptr == NULL)
        return -1;
    *lineptr = final_lineptr;
    *n = final_size;
    return (int) pos;
}

/* Handles errors by freeing all allocated memory and printing an error message.
 * Arguments:
 * n_mat, A, D, W : matrices to free
 * n              : number of rows in the matrices
 */
void handle_error_and_exit(double **n_mat, double **A, double **D, double **W, int n) {
    free_all(n_mat, A, D, W, n);
    printf("An Error Has Occurred\n");
    exit(1);
}

/* Frees the memory allocated for a matrix.
 * Arguments:
 * matrix : matrix to free
 * rows   : number of rows in the matrix
 */
void free_matrix(double **matrix, int rows) {
    int i;
    if (matrix == NULL) return;
    for (i = 0; i < rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

/* Frees all dynamically allocated matrices.
 * Arguments:
 * n_mat, A, D, W : matrices to free
 * n              : number of rows in the matrices
 */
void free_all(double **n_mat, double **A, double **D, double **W, int n) {
    free_matrix(n_mat, n);
    free_matrix(A, n);
    free_matrix(D, n);
    free_matrix(W, n);
}

/* Prints a matrix with each element formatted to four decimal places.
 * Arguments:
 * matrix : matrix to print
 * rows   : number of rows in the matrix
 * cols   : number of columns in the matrix
 */
void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols - 1; j++) {
            printf("%.4f,",
                   matrix[i][j]);  /* Print each element with format specifier for double separated by commas. */
        }
        printf("%.4f", matrix[i][j]); /* final number in each row doesn't need a comma */
        printf("\n");  /* New line after each row */
    }
}

/* Calculates the squared Euclidean distance between two vectors.
 * Arguments:
 * u, v : vectors of dimension d
 * d    : dimension of the vectors
 * Returns the squared Euclidean distance.
 */
double distance_squared(double *u, double *v, int d) {
    double res = 0;
    int i;
    for (i = 0; i < d; i++) {
        double deltapow = pow((u[i] - v[i]), 2);
        res += deltapow;
    }
    return res;
}

/* Calculates the symmetric matrix A.
 * Arguments:
 * n_mat : matrix of input vectors
 * n     : number of vectors
 * d     : dimension of vectors
 * flag  : pointer to error flag, set to 1 on failure
 * Returns the symmetric matrix.
 */
double** calc_sym_matrix(double **n_mat, int n, int d, int *flag) {

    int i, j;
    double curr_distance;
    double **A;
    A = (double **) calloc(n, sizeof(double *));
    if (A == NULL) {
        *flag = 1;
        return A;
    }
    for (i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
        if (A[i] == NULL) {
            *flag = 1;
            return A;
        }
        for (j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = 0;
            else {
                curr_distance = distance_squared(n_mat[i], n_mat[j], d);
                A[i][j] = exp(-curr_distance / 2);
            }
        }
    }
    return A;
}

/* Calculates the diagonal degree matrix D.
 * Arguments:
 * A     : symmetric matrix
 * n     : number of vectors
 * flag  : pointer to error flag, set to 1 on failure
 * Returns the diagonal matrix D.
 */
double** calc_ddg_matrix(double **A, int n,int *flag) {

    int i, j;
    double curr_sum;
    double **D;
    D = (double **) calloc(n, sizeof(double *));
    if (D == NULL) {
        *flag = 1;
        return D;
    }
    for (i = 0; i < n; i++) {
        D[i] = (double *) calloc(n, sizeof(double));
        if (D[i] == NULL) {
            *flag = 1;
            return D;
        }
        curr_sum = 0;
        for (j = 0; j < n; j++) {
            curr_sum += A[i][j];
        }
        D[i][i] = curr_sum;
    }
    return D;
}

/* Calculates the normalized matrix W.
 * Arguments:
 * A     : symmetric matrix
 * D     : diagonal degree matrix
 * n     : number of vectors
 * flag  : pointer to error flag, set to 1 on failure
 * Returns the normalized similarity matrix W.
 */
double** calc_norm_matrix(double **A, double **D, int n,int *flag) {

    int i, j;
    double **W;
    W = (double **) calloc(n, sizeof(double *));
    if (W == NULL) {
        *flag = 1;
        return W;
    }
    for (i = 0; i < n; i++) {
        W[i] = (double *) malloc(n * sizeof(double));
        if (W[i] == NULL) {
            *flag = 1;
            return W;
        }
        for (j = 0; j < n; j++) {
            W[i][j] = (1 / sqrt(D[i][i])) * A[i][j] * (1 / sqrt(D[j][j]));
        }
    }
    return W;
}
