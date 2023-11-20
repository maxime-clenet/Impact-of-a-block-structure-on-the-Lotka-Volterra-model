# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet

The aim is to display a dynamics of the LV system in
the case of a block matrix.
We can consider different communities and plot them
distinctly.
"""

# Importation of the important packages and main functions.
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
#from functions import block_matrix_normal, dynamics_LV_2B


def alpha_type_linear(t, N):
    if 15 > t > 5:
        return (t-5)*N
    elif t >= 15:
        return 100
    else:
        return 1
    

def alpha_type_normal(t, N):
    return 1



# Choix nombre itération:
nbr_it = 10000
# Choix du pas de temps de notre schéma:
tau = 0.002
# Choix condition initiale


N = 10

x_init = np.random.random(N)*2

x_init[:5] *= 1

Beta = [1/2, 1/2]
Sigma = np.array([[1/2, 1/80], [1/80, 1/2]])
Mu = np.array([[0, 0], [0, 0]])

A = block_matrix_normal(N, Beta, Mu, Sigma)
A_base = A.copy()

S = dynamics_LV_2B(alpha_type_normal, A=A, x_init=x_init,
                   N=N, nbr_it=nbr_it, tau=tau)


x = np.linspace(0, nbr_it*tau, nbr_it)

fig = plt.figure(1, figsize=(10, 6))

for i in range(S.shape[0]):
    #lab = '$x_'+str(i)+'$'
    if (i < S.shape[0]//2 and i == 0):
        plt.plot(x, S[i, :], label='Community 1', color='red')
    elif(i < S.shape[0]//2 and i != 0):
        plt.plot(x, S[i, :], color='red')
    elif(i == 5):
        plt.plot(x, S[i, :], label='Community 2', color='blue')
    else:
        plt.plot(x, S[i, :], color='blue')


plt.xlabel("Time (t)", fontsize=15)
plt.ylabel("Abundance ($x_i$)", fontsize=15)
plt.legend(loc='upper left', fontsize=15)
plt.show()
plt.close()



BMN = (A_base/np.sqrt(N)).copy()
t = 15
# BMN[5:, :5] = BMN[5:, :5]*alpha_type_linear(t, 10)
# BMN[:5, 5:] = BMN[:5, 5:]*alpha_type_linear(t, 10)

BMN[5:, :5] = BMN[5:, :5]*alpha_type_normal(t, 10)
BMN[:5, 5:] = BMN[:5, 5:]*alpha_type_normal(t, 10)

# plt.colorbar(plt.pcolor(BMN))

#sns.set_theme()

ax = sns.heatmap(BMN, vmin=-1.5, vmax=1.5, square=True, cmap="seismic",
                 xticklabels=np.arange(1, 11, 1), yticklabels=np.arange(1, 11, 1))

