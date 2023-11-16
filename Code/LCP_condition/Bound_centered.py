# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:21:35 2022

@author: Maxime
"""


import numpy as np
import matplotlib.pyplot as plt
from Sensibility_LCP.functions import block_matrix_normal


import numpy as np
#from create_block_matrix import block_matrix_normal


def block_matrix_bound(n_size, beta, mu, sigma):
    """
    Creation of the block matrix associated

    Parameters
    ----------
    n_size : int,
        dimension of the matrix.
    beta : list (matrix),
        correspond to the size of each block.
    mu : list (matrix),
        correspond to the mean of each block.
    sigma : list (matrix)
        correspond to the variance of each block.

    Returns
    -------
    A : list (matrix)
        The filled block matrix.

    N_SIZE = 1000
    beta = [1/2,1/2]
    sigma = np.array([[1,1/4],[1/4,1]])/np.sqrt(N_SIZE)

    mu = np.array([[0,4],[4,0]])/N_SIZE

    BMN = block_matrix_normal(1000,beta,mu,sigma)
    """

    B = np.size(beta)

    A = np.ones((
        int(n_size*beta[0]), int(n_size*beta[0])))*sigma[0, 0]+mu[0, 0]

    for i in range(1, B):
        A_bis = np.ones((
            int(n_size*beta[0]), int(n_size*beta[i])))*sigma[0, i]+mu[0, i]
        A = np.concatenate([A, A_bis], axis=1)

    for j in range(1, B):
        Aj = np.ones((
            int(n_size*beta[j]), int(n_size*beta[0])))*sigma[j, 0]+mu[j, 0]
        for k in range(1, B):
            A_bisj = np.ones((
                int(n_size*beta[j]), int(n_size*beta[k])))*sigma[j, k]+mu[j, k]
            Aj = np.concatenate([Aj, A_bisj], axis=1)

        A = np.concatenate([A, Aj], axis=0)

    return A


# %%


N = 100

beta = np.array([3/4, 1/4])
sig_1 = np.array([[1/2, 1/5], [1/3, 1/3]])
sigma_norm_1 = sig_1/np.sqrt(N)

mu = np.array([[0, 0], [0, 0]])/N


BMN_1 = block_matrix_normal(N, beta, mu, sigma_norm_1)

BMN_2 = block_matrix_bound(N, beta, mu, sig_1**2/N)

eig_BMN_1 = np.linalg.norm(BMN_1+BMN_1.transpose(), ord=2)

print(eig_BMN_1)

# print(2*np.linalg.norm(BMN_2+BMN_2.transpose(),ord=2)**(1/2))
#print(2*(np.linalg.norm((sig_1**2+sig_1.T**2)@np.diag(beta), ord=2))**(1/2))
print(2*(np.linalg.norm(np.diag(beta)**(1/2) @
      (sig_1**2+sig_1.T**2)@np.diag(beta)**(1/2)))**(1/2))
