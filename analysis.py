import sys
import numpy as np
from sklearn.metrics import silhouette_score
from symnmf import get_n_mat_final_H_tuple
from kmeans import get_kmeans_clusters


# Generates cluster labels for each data point based on k-means results.
# Parameters:
# - k_mat : list
#     List of clusters from k-means, where each cluster contains assigned vector indices
# - n : int
#     Total number of data points
# Returns:
# - np.array
#     Cluster labels for each data point, with label values representing cluster indices
def get_kmeans_labels(k_mat,n):
    labels = [-1] * n
    for index in range(len(k_mat)):
        cluster = k_mat[index]
        for assigned_vec_index in cluster.assigned_vectors_indexes:
            labels[assigned_vec_index] = index
    return np.array(labels)


# Generates cluster labels for each data point from SymNMF results
# Parameters:
# - H : np.array
#     Final matrix H from SymNMF, where each row represents a data point
# Returns:
# - np.array
#     Cluster labels for each data point, based on the index of the maximum value in each row.
def get_symnmf_labels(H):
    # Get the index of the maximum value in each row (axis=1)
    labels = np.argmax(H, axis=1)
    return labels

# Main function to compare clustering results between SymNMF and k-means using silhouette scores.
# Reads input parameters from the command line:
# - k : int
#     Number of clusters
# - file_name : str
#     File name  containing the data matrix
# For each clustering method, calculates labels and computes the silhouette score
# Prints silhouette scores for SymNMF and k-means
def main():
    eps = 0.0001
    maxIter = 300
    k = sys.argv[1]
    k = int(float(k))
    file_name = sys.argv[2]
    try:
        n_mat, final_H = get_n_mat_final_H_tuple(file_name, k)
        k_mat = get_kmeans_clusters(file_name, k, maxIter, eps)[1]
    except:
        print("An Error Has Occurred")
        return
    symnmf_labels = get_symnmf_labels(np.array(final_H))
    symnmf_score = silhouette_score(n_mat, symnmf_labels)
    print(f"nmf: {symnmf_score:.4f}")

    kmeans_labels = get_kmeans_labels(k_mat, len(n_mat))
    kmeans_score = silhouette_score(n_mat, kmeans_labels)
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()