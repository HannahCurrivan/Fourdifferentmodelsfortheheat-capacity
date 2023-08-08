import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

NA = 6.022e23 # Avogadros Number
BC = 1.38e-23 # Boltzmann Constant
PC = 1.055e-34 # Planck Constant
R = NA*BC # Gas Constant

# The four different models for the heat capacity

# Return heat capacity using Dulong Petit
def D_P(N, T):
    dulong_petit = 3*N*BC*np.ones(len(T))
    return dulong_petit

# Return heat capacity using Einstein Model
def einstein(N, ET, T):
    Einstein_model = 3*N*BC*(ET/T)**2*np.exp(ET/T)/(np.exp(ET/T)-1)**2
    return Einstein_model

# Return heat capacity using Debye Model
def debye(N, DT, T):
    C = np.zeros(len(T))
    for i, t in enumerate(T):
        Dx = DT/t
        I = quad(lambda x: x**4*np.exp(x)/(np.exp(x)-1)**2, 0, Dx)[0]
        C[i] = 9*N*BC*(t/DT)**3*I
        return C

# Return heat capacity using Debye Model Low Temperature Limit
def D_LowTemp(N, DT, T):
    Debye_LowTemp = 12*np.pi**4/5*N*BC*(T/DT)**3
    return Debye_LowTemp

# Specific quantities 
Tv = np.linspace(1, 500, 1000) # Array of temperature values
NO = NA # Total number of oscillators
Debye_T = 105 # Debye temperature
Einstein_T = Debye_T # Einstein temperature

# Plot the results
plt.figure()
plt.plot(Tv/Debye_T, D_P(NO,Tv), 'k-')
plt.plot(Tv/Debye_T, einstein(NO,Einstein_T,Tv), 'r')
plt.plot(Tv/Debye_T, debye(NO,Einstein_T,Tv), 'b')
plt.plot(Tv/Debye_T, D_LowTemp(NO, Einstein_T, Tv), 'b:')
plt.ylim([0,50])
plt.xlim([0, Tv[-1]/Debye_T])
plt.legend([r'Dulong-Petit', r'Einstein', r'Debye', r'Debye, low $T$'], loc=4)
plt.xlabel(r'$T/DT$')
plt.ylabel(r'$C_V \mathrm{[J/K]}$')
plt.show()
