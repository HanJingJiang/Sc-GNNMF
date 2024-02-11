import sys
import getopt
import numpy as np
import pandas as pd
from sklearn import metrics
import sklearn.cluster as cl


def POCR(inputfile, outputfile, cluster_num):
    # read matrix file as dataframe, the input file is preprocessed by R package
    input_df = pd.read_csv(inputfile, sep=",", header=None)
    n_cells, m_genes = input_df.shape
    input_df = np.log(input_df + 1)  # for every number x in matrix, x=log(x+1)
    if n_cells >= 400:
        print("big dataset. Use the linear kernel embedding cell-cell similarity")
        S = SLKE_similarity_linear(input_df, 1e-5, 0.001)
        labels = get_spectralClustering(S, cluster_num)
    else:
        print("small dataset. Use the RBF kernel embedding cell-cell similarity")
        S = SLKE_similarity_gauss(input_df, 1, 1e-5, 0.001)
        labels = get_spectralClustering(S, cluster_num)
    np.savetxt(outputfile, labels, fmt='%f', delimiter=',', encoding='utf-8')
    return True


def SLKES(K, a, b):
    n = K.shape[0]
    W = np.eye(n)
    J = W
    Z = W
    E = np.eye(n)
    Y1 = Y2 = np.zeros(n)
    count = 0
    D = K
    while (count < 600):
        Zold = Z.copy()
        J = np.linalg.inv(b * E + D * W * W.T * D.T) * (b * Z + Y1 + D * W * K.T)
        J[J < 0] = 0
        W = np.linalg.inv(b * E + D.T * J * J.T * D) * (b * Z + Y2 + D.T * J * K)
        W[W < 0] = 0

        H = (W - Y1 / b + J - Y2 / b) / 2

        for i in range(n):
            for j in range(n):
                Z[i, j] = max(abs(H[i, j]) - (a / (2 * b)), 0) * np.sign(H[i, j])

        Z[Z < 0] = 0

        Y1 = Y1 + b * (Z - J)
        Y2 = Y2 + b * (Z - W)

        count = count + 1
        b = b * 1.25

        p = np.linalg.norm(Zold, 'fro')
        q = np.linalg.norm(Z - Zold, 'fro')
        if (count > 5 and q < p * 1e-7):
            break

    return Z


def linear_kernel(X):
    kernel = np.dot(X, X.T)
    return kernel


def RBF_kernel(X, t):
    '''
    Gauss kernel
    :param X: n samples, m genes
    :return: gauss kernel
    '''
    n, m = X.shape
    X = pd.DataFrame(X)
    EU = metrics.euclidean_distances(X)
    d_max = EU.max()
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            result[i][j] = np.exp(-(EU[i][j] * EU[i][j]) / (t * d_max * d_max))
            result[j][i] = result[i][j]
    max = result.max()
    result = result / max
    return result


def SLKE_similarity_gauss(X, t, gamma, mu):
    K = RBF_kernel(X, t)
    K = matrix_normalize(K)
    S = SLKES(K, gamma, mu)
    S = matrix_normalize(S)
    return S


def SLKE_similarity_linear(X, gamma, mu):
    K = linear_kernel(X)
    K = matrix_normalize(K)
    S = SLKES(K, gamma, mu)

    S = matrix_normalize(S)
    return S


def matrix_normalize(similarity_matrix):
    similarity_matrix = np.matrix(similarity_matrix)
    similarity_matrix[np.isnan(similarity_matrix)] = 0
    if similarity_matrix.shape[0] == similarity_matrix.shape[1]:
        for i in range(similarity_matrix.shape[0]):
            similarity_matrix[i, i] = 1
        for i in range(200):
            D = np.diag(np.array(np.sum(similarity_matrix, axis=1)).flatten())
            D = np.linalg.pinv(np.sqrt(D))
            similarity_matrix = D * similarity_matrix * D
    else:
        for i in range(similarity_matrix.shape[0]):
            if np.sum(similarity_matrix[i], axis=1) == 0:
                similarity_matrix[i] = similarity_matrix[i]
            else:
                similarity_matrix[i] = similarity_matrix[i] / np.sum(similarity_matrix[i], axis=1)
    return similarity_matrix


def get_spectralClustering(similarity, cluster_num):
    """

    :param similarity: similarity matrix
    :param cluster_num: number of clusters(if it is 0, calculate by spectral clustering)
    :return: labels...
    """
    similarity = pd.DataFrame(similarity)
    similarity = similarity.values
    similarity[np.isnan(similarity)] = 0

    labels = cl.spectral_clustering(affinity=similarity, n_clusters=cluster_num)
    return labels


if __name__ == '__main__':
    inputfile = 'data/biase.csv'
    outputfile = 'labelbiase.csv'
    labelname = r'label//biase.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("biasefinished.")

    inputfile = 'data/yan.csv'
    outputfile = 'labelyan.csv'
    labelname = r'label//yan.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("yanfinished.")

    inputfile = 'data/usoskin.csv'
    outputfile = 'labelusoskin.csv'
    labelname = r'label//usoskin.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("usoskinfinished.")

    inputfile = 'data/pollen.csv'
    outputfile = 'labelpollen.csv'
    labelname = r'label//pollen.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("pollenfinished.")

    inputfile = 'data/li.csv'
    outputfile = 'labelli.csv'
    labelname = r'label//li.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("lifinished.")

    inputfile = 'data/10x-5cl.csv'
    outputfile = 'label10x-5cl.csv'
    labelname = r'label//10x-5cl.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("kolofinished.")

    inputfile = 'data/goolam.csv'
    outputfile = 'labelgoolam.csv'
    labelname = r'label//goolam.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("goolamfinished.")

    inputfile = 'data/deng.csv'
    outputfile = 'labeldeng.csv'
    labelname = r'label//deng.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("dengfinished.")

    inputfile = 'data/camp1.csv'
    outputfile = 'labelcamp1.csv'
    labelname = r'label//camp1.csv'
    with open(labelname, encoding='utf-8') as f:
        y = np.loadtxt(f, delimiter=",")
    y = y.astype(np.int)
    cluster_num = y.max() - y.min() + 1
    POCR(inputfile, outputfile, int(cluster_num))
    print("campfinished.")