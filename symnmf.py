import math
import sys
import numpy as np
import mysymnmfsp

# Set random seed for reproducibility
np.random.seed(1234)

# Define constants for maximum iterations, beta, and epsilon
maxIter = 300
beta = 0.5
eps = 0.0001

# Loads data from a text file and converts it to a list of lists (matrix format).
# Arguments:
# input_name : str - The filename of the input text file.
# Returns:
# List of lists with the data loaded from the file.
def init_n_mat(input_name):
    return np.loadtxt(input_name, delimiter=',').tolist() # returns data from input as np.array type


# Prints a matrix with each element formatted to four decimal places.
# Arguments:
# matrix : list of lists - The matrix to print.
def print_matrix(matrix):
    for v in matrix:
        res = ""
        for num in v:
            res = f"{res},{num:.4f}"  # in c maybe need to add \n
        print(res[1:])

# Calculates the mean value of a matrix.
# Arguments:
# W : list of lists - The input matrix.
# Returns:
# The mean value of the matrix elements as a float.
def calc_m(W):
    W = np.array(W)
    m = np.mean(W)
    return m

# Initializes the matrix H0 with random values within a specified range based on W.
# Arguments:
# W : list of lists - The input matrix used to determine the range.
# k : int - The number of columns for the initialized matrix H0.
# Returns:
# A list of lists representing matrix H0 with dimensions n x k.
def init_H(W, k):
    m=calc_m(W)

    # Define the matrix dimensions
    rows = len(W) # equals n
    cols = k

    # Define the interval [a, b]
    a = 0
    b = 2 * math.sqrt(m / k)

    # Generate a matrix with random values between a and b
    H0 = np.random.uniform(low=a, high=b, size=(rows, cols))
    return H0.tolist()

# Retrieves the normalized matrix W using the C extension function and computes the final matrix H using SymNMF.
# This function is used by analysis.py
# Arguments:
# file_name : str - The filename of the input CSV file.
# k : int - The number of clusters for SymNMF.
# Returns:
# Tuple containing the input matrix (n_mat) and the final computed matrix (H).
def get_n_mat_final_H_tuple(file_name,k):
    n_mat = init_n_mat(file_name)
    W = mysymnmfsp.norm(n_mat)
    H0 = init_H(W, k)
    H = mysymnmfsp.symnmf(H0, W, maxIter, beta, eps)
    return n_mat, H

# Main function to execute based on command-line arguments.
# The program expects three command-line arguments:
# - k : Number of clusters (int)
# - goal : Goal of the computation ("sym", "ddg", "norm", or "symnmf")
# - file_name : The filename of the input data
def main():
    k = sys.argv[1]
    k = int(float(k))
    goal = sys.argv[2]
    file_name = sys.argv[3]
    try:
        n_mat = init_n_mat(file_name)
        if goal == "sym":
            A = mysymnmfsp.sym(n_mat)
            print_matrix(A)

        elif goal == "ddg":
            D = mysymnmfsp.ddg(n_mat)
            print_matrix(D)

        elif goal == "norm":
            W = mysymnmfsp.norm(n_mat)
            print_matrix(W)

        elif goal == "symnmf":
            H = get_n_mat_final_H_tuple(file_name, k)[1]
            print_matrix(H)

    except:
        print("An Error Has Occurred")

if __name__ == "__main__":
    main()
