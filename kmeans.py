import sys
import math


# Class for a single cluster
class Cluster:

    # Initializes the cluster with a given vector as its initial centroid.
    # Arguments:
    # vector: List of floats, representing the initial centroid for the cluster.
    def __init__(self, vector):
        self.centroid = [num for num in vector]
        self.prev_centroid = None  # Holds the previous centroid for convergence checking
        self.d = len(vector)  # Dimension of the vector
        self.sum_vector = [0] * self.d  # Sum of vectors assigned to this cluster
        self.size = 0  # Number of vectors assigned to this cluster
        self.assigned_vectors_indexes = []  # Indexes of vectors assigned to this cluster

    # Resets assigned vectors, size, and sum vector for a new iteration
    def reset(self):
        self.sum_vector = [0] * self.d
        self.size = 0
        self.assigned_vectors_indexes = []
    # Adds a vector to the cluster by updating the sum and increasing the size
    # Arguments:
    # x_vector: List of floats, representing the vector to add.
    # index: Integer, index of the vector in the main dataset.
    def add(self, x_vector,index):
        self.size += 1
        self.assigned_vectors_indexes.append(index)
        for i in range(self.d):
            self.sum_vector[i] += x_vector[i]


    # Updates the cluster's centroid based on the average of assigned vectors.
    def update(self):
        self.prev_centroid = [num for num in self.centroid]
        self.centroid = [(sum / self.size) for sum in self.sum_vector]


    # Returns detailed information about the cluster, including centroid and assigned vectors.
    # Arguments:
    # n_mat: List of lists, where each list is a vector in the dataset.
    def get_full_info(self,n_mat):
        res="Cluster: " + str(self.centroid) +  "\n"
        res+="Assigned Vectors:"
        for assigned_vec_index in self.assigned_vectors_indexes:
            res += "\nVector index = " + str(assigned_vec_index)
            res += ", Vector = " + str(n_mat[assigned_vec_index])
        return res


    # Returns a string representation of the cluster's centroid with values formatted to four decimal places.
    def __repr__(self):
        res = ""
        for num in self.centroid:
            res = f"{res},{num:.4f}"    #in c maybe need to add \n
        return res[1:]


# Calculates the Euclidean distance between two vectors.
# Arguments:
# u, v: Lists of floats representing vectors of the same dimension.
def distance(u, v):
    n = len(u)
    if n != len(v):
        print("An Error Has Occurred")
        return -1

    sum_under_root = 0
    for i in range(n):
        sum_under_root += (u[i] - v[i]) ** 2
    return math.sqrt(sum_under_root)


# Checks if a number is a valid integer.
# Arguments:
# num: The number to check.
# Returns 0 if the number is an integer, 1 otherwise.
def is_int(num):
    try:
        num=float(num)
        if num != int(num):
            return 1
        return 0
    except:
        return 1


# Validates numeric inputs for k and iter.
# Arguments:
# k, iter: The values to check for validity.
# Returns False if either value is invalid, True otherwise.
def check_numeric_validity(k, iter):
    if is_int(k) == 1:
        print("Invalid number of clusters!")
        return False
    if is_int(iter) == 1:
        print("Invalid maximum iteration!")
        return False
    return True


# Validates that k and iter are within acceptable ranges.
# Arguments:
# k, iter: The number of clusters and the max iterations.
# N: Total number of vectors in the dataset.
def check_range_validity(k, iter, N):
    if k <= 1 or k >= N:
        print("Invalid number of clusters!")
        return False
    if iter <= 1 or iter >= 1000:
        print("Invalid maximum iteration!")
        return False
    return True


# Reads lines from a file and returns them as a list of strings.
# Arguments:
# input: The file path to read from.
def get_lines_from_file(input):
    try:
        f = open(input, "r")
        input_lines = [line.strip() for line in f.readlines()]
        f.close()
        return input_lines
    except:
        print("An Error Has Occurred")
        return -1


# Initializes the dataset matrix (n_mat) from input lines.
# Arguments:
# input_lines: List of strings representing lines of comma-separated floats.
def init_n_mat(input_lines):
    res = []
    for line in input_lines:
        res.append([float(num) for num in line.split(",")])
    return res


# Initializes the clusters with the first k vectors as initial centroids.
# Arguments:
# n_mat: List of vectors in the dataset.
# k: Number of clusters to create.
def init_k_mat(n_mat, k):
    k_mat = []
    for i in range(k):
        k_mat.append(Cluster(n_mat[i]))
    return k_mat


# Finds the cluster with the minimum distance to a given vector.
# Arguments:
# k_mat: List of clusters.
# x_vector: Vector to find the closest cluster for.
def find_min_cluster(k_mat, x_vector):
    min_cluster = k_mat[0]
    min_d = distance(x_vector, k_mat[0].centroid)  # maybe k_mat is empty?
    for i in range(1, len(k_mat)):
        curr_d = distance(x_vector, k_mat[i].centroid)
        if curr_d < min_d:
            min_cluster = k_mat[i]
            min_d = curr_d
    return min_cluster


# Updates centroids of all clusters based on current assigned vectors.
# Arguments:
# k_mat: List of clusters.
def update_clusters(k_mat):
    for cluster in k_mat:
        cluster.update()


# Resets all clusters to prepare for a new iteration.
# Arguments:
# k_mat: List of clusters.
def reset_clusters(k_mat):
    for cluster in k_mat:
        cluster.reset()


# Checks if all clusters have converged based on a specified epsilon.
# Arguments:
# k_mat: List of clusters.
# epsilon: Threshold for convergence.
def check_conv(k_mat,epsilon):
    for cluster in k_mat:
        if distance(cluster.centroid, cluster.prev_centroid) >= epsilon:
            return False
    return True


# Runs the k-means algorithm and returns the dataset and clusters.
# Arguments:
# file_name: Path to input file with vectors.
# k: Number of clusters.
# maxIter: Maximum iterations for the algorithm.
# epsilon: Convergence threshold.
def get_kmeans_clusters(file_name,k,maxIter,epsilon):
    conv = False
    input_lines = get_lines_from_file(file_name)
    if input_lines == -1:
        return None
    N = len(input_lines)
    if not check_range_validity(k, maxIter, N):
        return None
    n_mat = init_n_mat(input_lines)
    k_mat = init_k_mat(n_mat, k)
    iter_count = 1
    while iter_count <= maxIter and not conv:
        reset_clusters(k_mat)
        for index in range(len(n_mat)):
            x_vector=n_mat[index]
            min_cluster = find_min_cluster(k_mat, x_vector)
            # presumes k_mat has been reset.
            min_cluster.add(x_vector,index)
        update_clusters(k_mat)
        iter_count += 1
        conv = check_conv(k_mat,epsilon)
    return n_mat, k_mat


# Main function to run the k-means algorithm with command line arguments.
def main():
    eps=0.001
    k = sys.argv[1]
    file_name = sys.argv[-1]
    maxIter = sys.argv[2] if len(sys.argv) == 4 else 200
    if not check_numeric_validity(k, maxIter):
        return 1
    k = int(float(k))
    maxIter = int(float(maxIter))
    n_mat,k_mat=get_kmeans_clusters(file_name,k,maxIter,eps)
    for cluster in k_mat:
        print(cluster.get_full_info(n_mat))


if __name__ == "__main__":
    main()
