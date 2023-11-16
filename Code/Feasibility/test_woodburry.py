# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 19:06:44 2022

@author: maxim
"""


import numpy as np
import matplotlib.pyplot as plt
from functions import block_matrix_normal

N = 10

beta = np.array([1/2, 1/2])
sig_1 = np.array([[0, 0], [0, 0]])
sigma_norm_1 = sig_1/np.sqrt(N)

mu = np.array([[1, 1], [1, 1]])

mu_norm = mu/N

I = np.eye(N)
BMN_1 = block_matrix_normal(N, beta, mu_norm, sigma_norm_1)

eig_BMN_1 = np.linalg.eig(I-BMN_1)[1]
eig = np.linalg.eig(BMN_1)[0]

inv = np.linalg.inv(I-BMN_1)

np.linalg.inv(np.linalg.inv(mu)-np.diag(beta))

np.linalg.inv((np.eye(2)-np.diag(beta)@mu@np.diag(beta)))


np.linalg.inv(np.linalg.inv(mu)-np.diag(beta))/10

A = -np.diag(beta)

(np.linalg.inv(A)-np.linalg.inv(A+A@mu@A))/10

(mu-mu@np.linalg.inv(np.linalg.inv(np.diag(-beta))+mu))/10
