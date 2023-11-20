# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet
"""

import numpy as np
import matplotlib.pyplot as plt


def f(kappa_12, kappa_11=1.2):

    a = 1-2/kappa_12**2
    b = 2/kappa_11**2-2/kappa_12**2
    if a/b > 0 and a/b < 1:
        return a/b
    elif a/b > 1:
        return 1
    else:
        return 0

def g(kappa_21, kappa_22=1.4):
    a = 1-2/kappa_22**2
    b = 2/kappa_21**2-2/kappa_22**2
    if a/b > 0 and a/b < 1:
        return a/b
    elif a/b > 1:
        return 1
    else:
        return 0



x = np.linspace(np.sqrt(2), 3, 100)
y_1 = np.linspace(0, 1, 100)
y_2 = np.linspace(0, 1, 100)

for i in range(len(x)):
    y_1[i] = f(x[i])
    y_2[i] = g(x[i])

fig = plt.figure(1, figsize=(10, 6))
plt.plot(x, y_1, label='Upper bound', color="#0000ff")
plt.plot(x, y_2, label='Lower bound', color='red')
plt.fill_between(x, y_1, y_2, where=y_1 > y_2, color='#4300c84d')
plt.xlabel(r"$\kappa_{12},\kappa_{21}$", fontsize=15)
plt.ylabel(r"$\beta_1$", fontsize=15)
plt.legend(loc='upper right')
plt.show()

