# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet


This file allows you to display the graph showing the
theoretical distribution vs. the empirical distribution.
"""

# Importation of the packages and required functions.

import numpy as np
import matplotlib.pyplot as plt
from lemkelcp import lemkelcp


# For more informations on this two functions see functions.py
from functions import *


def plot_distrib(mu, s, beta, B_size=2000, law_type='normal'):
    """
    Return a  graph showing the difference
    between theoretical distribution vs. the empirical distribution
    of the abondances.

    Parameters
    ----------
    alpha : float in [sqrt(2),Infty), optional
        Associated to the alpha value in the paper.
        The default is 2.
    mu : float in (-infty,1], optional
        Associated to the mu value in the paper.
        The default is 0.2.
    B_size : int, optional
        Size of the square matrix B.
        The default is 2000.

    Returns
    -------
    fig : matplotlib.figure.

    """
    # Computations for the empirical distribution.
    # We find a solution using the pivot algorithm.
    if law_type == 'normal':
        mu_norm = mu/B_size
        Sig = s
        sigma_norm = Sig/np.sqrt(B_size)
        B = block_matrix_normal(B_size, beta, mu_norm, sigma_norm)
    if law_type == 'unif':
        Sig = s
        sigma_norm = Sig/np.sqrt(B_size)
        B = block_matrix_unif(B_size, beta, sigma_norm)

    q = np.ones(B_size)
    M = -np.eye(B_size)+B

    res_lcp = lemkelcp.lemkelcp(-M, -q, maxIter=10000)[0]
    res_lcp_1 = res_lcp[:int(beta[0]*B_size)]
    res_lcp_2 = res_lcp[int(beta[0]*B_size):]
    res_lcp_pos_1 = res_lcp_1[res_lcp_1 != 0]
    res_lcp_pos_2 = res_lcp_2[res_lcp_2 != 0]
    (pi_1, pi_2, m_1, m_2, sigma_1, sigma_2) = block_function(mu, 1/s, beta)
    print(m_1,m_2)
    x = np.linspace(0.01, 4, 1000)

    fig = plt.figure(1, figsize=(10, 6))
    y_1 = np.ones(len(x))
    y_2 = np.ones(len(x))
    for i, v in enumerate(x):
        y_1[i] = block_density(v, 0, pi_1, pi_2, m_1, m_2,
                               sigma_1, sigma_2, mu, Sig, beta)
    
    for i, v in enumerate(x):
        y_2[i] = block_density(v, 1, pi_1, pi_2, m_1, m_2,
                               sigma_1, sigma_2, mu, Sig, beta)
    if law_type == 'normal':
        plt.vlines(np.mean(res_lcp_pos_1),0,3,linewidth=2,linestyles = "dashed",color = "blue")
        plt.vlines(np.mean(res_lcp_pos_2),0,3,linewidth=2,linestyles = "dashed",color = "red")
        plt.ylim(0,2.3)

    plt.plot(x, y_1, linewidth=2.5, color='blue', label='Community 1')
    plt.plot(x, y_2, linewidth=2.5, color='red', label='Community 2')
    plt.hist(res_lcp_pos_1, density=True, bins=20,
             edgecolor='blue', color='#0000003d')
    plt.hist(res_lcp_pos_2, density=True, bins=20,
             edgecolor='red', color='#0000003d')
    plt.xlabel('Abundances (x*)', fontsize=15)
    plt.ylabel('Density of the distribution (f)', fontsize=15)
    plt.xlim(0, 4)

    plt.legend(loc='upper right')
    plt.show()

    return fig


# mu = np.array([[0, 0], [0, 0]])
# s = np.array([[1/2, 1/np.sqrt(2)], [1/5, 1/9]])
# beta = np.array([3/4, 1/4])

# plot_distrib(mu, s, beta, B_size=2000, law_type='normal')


mu = np.array([[0, 0], [0, 0]])
s = np.array([[1/2, 2/3], [1/3, 1/4]])
beta = np.array([1/2, 1/2])

plot_distrib(mu, s, beta, B_size=2000, law_type='unif')

