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


![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/EM.JPG)

### Debye Model:


![alt text](https://github.com/HannahCurrivan/Fourdifferentmodelsfortheheat-capacity/blob/main/DM.JPG)

### Debye Model Low-Temperature Limit:


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
# Return heat capacity using Debye Model Low Temperature Limit
def D_LowTemp(N, DT, T):
    Debye_LowTemp = 12*np.pi**4/5*N*BC*(T/DT)**3
    return Debye_LowTemp
```
