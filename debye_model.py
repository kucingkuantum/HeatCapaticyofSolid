""" program to calculate 
    Heat capacity of solid
    by using Debye Model, Einstein
    and Dulong-Petit
    Written by Syahril Siregar
    2018/1/16
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def integrand(x):
    return (x**4)*np.exp(x)/((np.exp(x)-1)**2)

n = 1.0 # mole
R = 8.314 #J/K
Tetad = 215.0 #K
Tetae = 215.0
x0 = 1.0e-12 # 0 boundary
Tmin = 1   #K
Tmax = 600 #K
cvd = np.zeros(Tmax-Tmin)
cve = np.zeros(Tmax-Tmin)
cdp = np.zeros(Tmax-Tmin)
Temp = np.linspace(Tmin,Tmax, Tmax-Tmin, endpoint=True)

#Debye Calculation 
for T in range(Tmin,Tmax):
    a,b = quad(integrand, x0, Tetad/T,args=())
    cvd[T-Tmin] = a*n*R*9.0*((T/Tetad)**3)
#Einstein
    cve[T-Tmin] = 3*n*R*((Tetae/T)**2)*(np.exp(Tetae/T)/(np.exp(Tetae/T)-1)**2)
#Dulong-Petit
    cdp[T-Tmin] = 3*n*R

#plotting
l1=plt.plot(Temp,cvd)
l2=plt.plot(Temp,cve)
l3=plt.plot(Temp,cdp,'--')
plt.xlabel('Temperature(K)')
plt.ylabel('Heat Capacity (J/K mole) ')
plt.legend(('Debye Model', 'Einstein Model', 'Dulong-Petit'),
           loc='lower right')
plt.show()
