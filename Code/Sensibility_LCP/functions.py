# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 11:41:17 2022

@author: Maxime

"""


# Importation of the packages:

import numpy as np
from scipy import optimize
import scipy.stats as stats
from lemkelcp import lemkelcp


# %%

# The main functions used in the block model:

def block_matrix_normal(n_size, beta, mu, sigma):
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

    A = np.random.randn(
        int(n_size*beta[0]), int(n_size*beta[0]))*sigma[0, 0]+mu[0, 0]

    for i in range(1, B):
        A_bis = np.random.randn(
            int(n_size*beta[0]), int(n_size*beta[i]))*sigma[0, i]+mu[0, i]
        A = np.concatenate([A, A_bis], axis=1)

    for j in range(1, B):
        Aj = np.random.randn(
            int(n_size*beta[j]), int(n_size*beta[0]))*sigma[j, 0]+mu[j, 0]
        for k in range(1, B):
            A_bisj = np.random.randn(
                int(n_size*beta[j]), int(n_size*beta[k]))*sigma[j, k]+mu[j, k]
            Aj = np.concatenate([Aj, A_bisj], axis=1)

        A = np.concatenate([A, Aj], axis=0)

    return A


def block_matrix_unif(n_size, beta, sigma):

    B = np.size(beta)

    A = (np.random.random((
        int(n_size*beta[0]), int(n_size*beta[0])))*2 *
        np.sqrt(3)-np.sqrt(3))*sigma[0, 0]

    for i in range(1, B):
        A_bis = (np.random.random((
            int(n_size*beta[0]), int(n_size*beta[i])))*2 *
            np.sqrt(3)-np.sqrt(3))*sigma[0, i]
        A = np.concatenate([A, A_bis], axis=1)

    for j in range(1, B):
        Aj = (np.random.random((
            int(n_size*beta[j]), int(n_size*beta[0])))*2 *
            np.sqrt(3)-np.sqrt(3))*sigma[j, 0]
        for k in range(1, B):
            A_bisj = (np.random.random((
                int(n_size*beta[j]), int(n_size*beta[k])))*2 *
                np.sqrt(3)-np.sqrt(3))*sigma[j, k]
            Aj = np.concatenate([Aj, A_bisj], axis=1)

        A = np.concatenate([A, Aj], axis=0)

    return A

# %%

# Part dedicated to the dynamics of the 2-blocks model


def f_LV_2B(x, N, A, alpha, A_1, A_2):
    A[0:5, 5:10] = A_1*alpha
    A[5:10, 0:5] = A_2*alpha
    x = np.dot(np.diag(x), (np.ones(N)-x+np.dot(A*(1/(np.sqrt(N))), x)))

    return(x)


def dynamics_LV_2B(alpha, A, x_init, N, nbr_it, tau):

    x = x_init

    A_1 = A[0:5, 5:10].copy()
    A_2 = A[5:10, 0:5].copy()
    compt = 0

    # Matrix of the solution:
    S = np.eye(N, nbr_it)

    # Euler scheme: (explicite version)
    while (compt < nbr_it):

        t = tau*compt

        f1 = f_LV_2B(x, N, A, alpha(t, N), A_1, A_2)
        f2 = f_LV_2B(x+tau*0.5*f1, N, A, alpha(t, N), A_1, A_2)
        f3 = f_LV_2B(x+tau*0.5*f2, N, A, alpha(t, N), A_1, A_2)
        f4 = f_LV_2B(x+tau*f3, N, A, alpha(t, N), A_1, A_2)

        x = x+tau*(f1+2*f2+2*f3+f4)/6

        for i in range(N):
            S[i, compt] = x[i]
        compt = compt+1

    return S


# %%


def E_cond(delta):

    a = -delta

    p_1 = np.exp(-a**2/2)
    p_2 = 1-stats.norm.cdf(-a)

    return (1/np.sqrt(2*np.pi))*p_1/p_2


def E_2_cond(delta):
    a = -delta
    p_1 = np.exp(-a**2/2)
    p_2 = 1-stats.norm.cdf(-a)

    return (1/np.sqrt(2*np.pi))*-a*p_1/p_2+1


def h_sigma_1(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D1 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[0, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[0, 1]**2)
    Lamb1 = pi_1*m_1*Beta[0]*mu[0, 0]+pi_2*m_2*Beta[1]*mu[0, 1]
    delta_1 = (-1-Lamb1)/D1

    a1 = 1+Lamb1
    b1 = D1

    return a1+b1*E_cond(delta_1)-m_1


def h_sigma_2(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D2 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[1, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[1, 1]**2)
    Lamb2 = pi_1*m_1*Beta[0]*mu[1, 0]+pi_2*m_2*Beta[1]*mu[1, 1]
    delta_2 = (-1-Lamb2)/D2

    a2 = 1+Lamb2
    b2 = D2

    return a2+b2*E_cond(delta_2)-m_2


def g_sigma_1(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D1 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[0, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[0, 1]**2)
    Lamb1 = pi_1*m_1*Beta[0]*mu[0, 0]+pi_2*m_2*Beta[1]*mu[0, 1]
    delta_1 = (-1-Lamb1)/D1

    a1 = (1+Lamb1)**2
    b1 = 2*(1+Lamb1)*D1
    c1 = D1**2

    return a1+b1*E_cond(delta_1)+c1*E_2_cond(delta_1)-sigma_1**2


def g_sigma_2(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D2 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[1, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[1, 1]**2)
    Lamb2 = pi_1*m_1*Beta[0]*mu[1, 0]+pi_2*m_2*Beta[1]*mu[1, 1]
    delta_2 = (-1-Lamb2)/D2

    a2 = (1+Lamb2)**2
    b2 = 2*(1+Lamb2)*D2
    c2 = D2**2

    return a2+b2*E_cond(delta_2)+c2*E_2_cond(delta_2)-sigma_2**2


def f_sigma_1(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D1 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[0, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[0, 1]**2)
    Lamb1 = pi_1*m_1*Beta[0]*mu[0, 0]+pi_2*m_2*Beta[1]*mu[0, 1]
    delta_1 = (-1-Lamb1)/D1

    return 1-stats.norm.cdf(delta_1)-pi_1


def f_sigma_2(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    D2 = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[1, 0]
                 ** 2+pi_2*sigma_2**2*Beta[1]*Sig[1, 1]**2)
    Lamb2 = pi_1*m_1*Beta[0]*mu[1, 0]+pi_2*m_2*Beta[1]*mu[1, 1]
    delta_2 = (-1-Lamb2)/D2

    return 1-stats.norm.cdf(delta_2)-pi_2


def Gamma(x, mu, Sig, Beta):
    """


    Parameters
    ----------
    x : x[O] correspond à sigma
        x[1] correspond à pi
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    return (f_sigma_1(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta), f_sigma_2(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta), h_sigma_1(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta), h_sigma_2(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta), g_sigma_1(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta), g_sigma_2(x[0], x[1], x[2], x[3], x[4], x[5], mu, Sig, Beta))


# %%


def block_function(mu, Alpha, Beta):
    """


    Parameters
    ----------
    mu : TYPE
        DESCRIPTION.
    Sig : TYPE
        DESCRIPTION.
    Beta : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    Sig = 1/Alpha

    (pi_1, pi_2, m_1, m_2, sigma_1, sigma_2) = optimize.root(
        Gamma, [0.999, 0.9999, 2, 2, 1.001, 1.001], args=(mu, Sig, Beta,)).x

    return(pi_1, pi_2, m_1, m_2, sigma_1, sigma_2)

# %%

# Empirical resolution:


def zero_LCP(A, beta):
    """
    This function resolve the LCP problem of our model.
    If a solution exist, this function return the properties
    of the solution i.e:
    - proportion of persistent species,
    - variance of the persistent species,
    - mean of the persistent species.

    Parameters
    ----------
    A : numpy.ndarray(n,n),
        Corresponds to the matrix of interactions.

    Returns
    -------
    If this solution exists, the function return: (En pourcentage)
    I/N : Proportion of surviving species.

    m : Mean of the surviving species.

    sigma: Root mean square of the surviving species.

    """
    A_SIZE = A.shape[0]
    q = np.ones(A_SIZE)
    M = -np.eye(A_SIZE)+A
    sol = lemkelcp.lemkelcp(-M, -q, maxIter=10000)

    I1 = int(beta[0]*A_SIZE)
    I2 = A_SIZE-I1

    res_LCP = sol[0]
    res_LCP_1 = res_LCP[:I1]
    res_LCP_2 = res_LCP[I1:]

    res_LCP_pos_1 = res_LCP_1[res_LCP_1 != 0]
    res_LCP_pos_2 = res_LCP_2[res_LCP_2 != 0]

    S1 = len(res_LCP_pos_1)
    S2 = len(res_LCP_pos_2)

    m1 = sum(res_LCP_pos_1)/S1
    m2 = sum(res_LCP_pos_2)/S2

    sigma1 = np.sqrt(sum(res_LCP_pos_1**2)/S1)
    sigma2 = np.sqrt(sum(res_LCP_pos_2**2)/S2)

    return (S1/I1, S2/I2, m1, m2, sigma1, sigma2)


def empirical_prop(A_size, alpha, mu, beta, mc_prec=500):
    """
    For a large number of matrix (mc_prec) of size (B_size), an empirical
    estimator of the parameter are given using a MC experiment.

    """

    S_p1 = np.zeros(mc_prec)
    S_sigma1 = np.zeros(mc_prec)
    S_m1 = np.zeros(mc_prec)
    S_p2 = np.zeros(mc_prec)
    S_sigma2 = np.zeros(mc_prec)
    S_m2 = np.zeros(mc_prec)

    Sig = 1/alpha
    Sig_norm = Sig/np.sqrt(A_size)
    mu_norm = mu/A_size
    for j in range(mc_prec):

        A = block_matrix_normal(A_size, beta, mu_norm, Sig_norm)

        (S_p1[j], S_p2[j], S_m1[j], S_m2[j],
         S_sigma1[j], S_sigma2[j]) = zero_LCP(A, beta)

    return np.mean(S_p1), np.mean(S_p2), np.mean(S_m1), np.mean(S_m2), np.mean(S_sigma1), np.mean(S_sigma2)


# %%

def block_density(x, k, pi_1, pi_2, m_1, m_2, sigma_1, sigma_2, mu, Sig, Beta):

    Delta = np.sqrt(pi_1*sigma_1**2*Beta[0]*Sig[k, 0]
                    ** 2+pi_2*sigma_2**2*Beta[1]*Sig[k, 1]**2)
    Lamb = pi_1*m_1*Beta[0]*mu[k, 0]+pi_2*m_2*Beta[1]*mu[k, 1]
    delta = (-1-Lamb)/Delta

    a = np.exp(-(x/Delta+delta)**2/2)/Delta
    b = stats.norm.cdf(-delta)*np.sqrt(2*np.pi)
    return a/b
