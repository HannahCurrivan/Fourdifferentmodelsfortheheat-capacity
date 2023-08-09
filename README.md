# Four different models for the heat capacity

This code displays the four possible outcomes of heat capacity for Silicon (Si). 

These models included:

- Dulong Petit
- Einstein Model
- Debye Model
- Debye Model Low-Temperature Limit

## Thermodynamic Laws:
### Dulong Petit:

Taken from Britannica; "Dulong–Petit law, statement that the gram-atomic heat capacity (specific heat times atomic weight) of an element is a constant; that is, it is the same for all solid elements, about six calories per gram atom."

![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/DPL.JPG)

### Einstein Model:
Taken from chem.libretext; "Einstein assumed three things when he investigated the heat capacity of solids. First, he assumed that each solid was composed of a lattice structure consisting of  N atoms. Each atom was treated as moving independently in three dimensions within the lattice (3 degrees of freedom). This meant that the entire lattice's vibrational motion could be described by a total of  3N motions, or degrees of freedom. Secondly, he assumed that the atoms inside the solid lattice did not interact with each other and thirdly, all of the atoms inside the solid vibrated at the same frequency. The third point highlights the main difference in Einstein's and Debye's two models."

![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/EM.JPG)

### Debye Model:
Taken from LibreTexts; Debye model estimating the phonon contribution to the specific heat (heat capacity) in a solid. This model correctly explains the low-temperature dependence of the heat capacity, which is proportional to  T3
  and also recovers the Dulong-Petit law at high temperatures. 

![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/DM.JPG)

### Debye Model Low-Temperature Limit:

This model correctly explains the low-temperature dependence of the heat capacity, which is proportional to T3 and also recovers the Dulong-Petit law at high temperatures.

![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/DLTM.JPG)

## Using Python to calculate the heat capacity of silicon:

Python packages needed to run this script:

```
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
```

First, you need to list the constants used across the four models. In this case, it is Avogadro's number, Boltzmann Constant, Planck Constant, and Gas Constant.

```
NA = 6.022e23 # Avogadros Number
BC = 1.38e-23 # Boltzmann Constant
PC = 1.055e-34 # Planck Constant
R = NA*BC # Gas Constant
```

Now to apply the four different models, starting with Dulong Petit:

```
# Return heat capacity using Dulong Petit
def D_P(N, T):
    dulong_petit = 3*N*BC*np.ones(len(T))
    return dulong_petit
```

Applying the Einstein model to find Silicon heat capacity:

```
# Return heat capacity using Einstein Model
def einstein(N, ET, T):
    Einstein_model = 3*N*BC*(ET/T)**2*np.exp(ET/T)/(np.exp(ET/T)-1)**2
    return Einstein_model
```

Applying the Debye Model to find Silicon heat capacity:

```
# Return heat capacity using Debye Model
def debye(N, DT, T):
    C = np.zeros(len(T))
    for i, t in enumerate(T):
        Dx = DT/t 
        I = quad(lambda x: x**4*np.exp(x)/(np.exp(x)-1)**2, 0, Dx)[0]
        C[i] = 9*N*BC*(t/DT)**3*I
        return C
```

Applying the Debye Model Low-Temperature Limit to find Silicon heat capacity:

```
# Return heat capacity using Debye Model Low-Temperature Limit
def D_LowTemp(N, DT, T):
    Debye_LowTemp = 12*np.pi**4/5*N*BC*(T/DT)**3
    return Debye_LowTemp
```

 Now you have to apply the following specific quantities to each of the function's inputs.

```
# Specific quantities 
Tv = np.linspace(1, 500, 1000) # Array of temperature values
NO = NA # Total number of oscillators
Debye_T = 105 # Debye temperature
Einstein_T = Debye_T # Einstein temperature
```

Now you can plot the four different heat capacity models for Silicon. 

```
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
```

The graphical outcome:

![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/LDT.png)
