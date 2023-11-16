import numpy as np
import matplotlib.pyplot as plt

x = -np.linspace(-1,2,1000)
z = np.linspace(-2,2,1000)


beta1 = 1/2
beta2 = 1-beta1
gamma_11 = 0.5
gamma_22 = 0.5
Gamma_n = 0.8

c = Gamma_n-beta1*gamma_11**2-beta2*gamma_22**2

x,z = np.meshgrid(x,z)

u = 0  # x-position of the center
v = 0  # y-position of the center
a = np.sqrt(c/beta2)  # radius on the x-axis
print(a)
b = np.sqrt(c/beta1) # radius on the y-axis
t = np.linspace(0, 2*np.pi, 100)



fig = plt.figure(1, figsize=(9, 6))

plt.xlabel(r"$\gamma_{12}$", fontsize=15)
plt.ylabel(r"$\gamma_{21}$", fontsize=15)
plt.vlines(np.sqrt((1-2*beta1*gamma_11**2)/(2*beta2)),-2,4,linestyles = 'dashed',color='black')
plt.hlines(np.sqrt((1-2*beta2*gamma_22**2)/(2*beta1)),-2,4,linestyles = 'dotted',color='black')

#plt.hlines(np.sqrt((kappa-beta2*s22**2-1/2)/beta1),-2,4,linestyles = 'dashed')

plt.plot(u+a*np.cos(t), v+b*np.sin(t), color='black', linewidth=2)
plt.grid(color='lightgray', linestyle='--')
plt.axis("equal")
plt.xlim(0, 2.2)
plt.ylim(0, 1.5)
plt.xticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]) 
plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4]) 
plt.show()





