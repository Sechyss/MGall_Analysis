#%% Imports

import numpy as np
import pandas as pd
from scipy.datasets import face
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#%% Parameters

# Initial population values
S = 10000  # Initial susceptible population
Eh = 0  # Initial exposed humans high virulence pathogen
Indh = 5  # Initial infected humans high virulence pathogen not taking the drug
Idh = 0  # Initial infected humans high virulence pathogen taking the drug
Rh = 0  # Initial recovered humans high virulence pathogen
El = 0  # Initial exposed humans low virulence pathogen
Indl = 5  # Initial infected humans low virulence pathogen not taking the drug
Idl = 0  # Initial infected humans low virulence pathogen taking the drug
Rl = 0  # Initial recovered humans low virulence pathogen

pop_values = (S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl)
N = sum(pop_values)  # Total initial population


# Model parameters

theta = 0.6  # Proportion of population that takes the drug
p_mortality = 0.1  # Protection effect of the drug against mortality
p_recover = 1  # Protection effect of the drug for recovery
phi_transmission = 1.5  # Increase in transmission rate for high virulence pathogen
phi_mortality = 2  # Increase in mortality rate for high virulence pathogen
phi_recover = 1  # Decrease in recovery rate for high virulence pathogen
zeta = 0.02  # Death rate due to infection
sigma = 1/7  # Rate of recovery from infection
delta = 1/150  # Rate of loss of immunity
tau = 1/5  # Rate of progression from exposed to infected

death_rate = 0.01 # Natural death rate due to other causes
birth_rate = zeta + death_rate # Birth rate


beta_l = 0.3  # Transmission rate low virulence pathogen
beta_h = phi_transmission * beta_l  # Transmission rate high virulence pathogen

beta_values = beta_l, beta_h
parameters = (birth_rate, death_rate, theta, p_mortality, p_recover, tau,
              phi_mortality, phi_recover, zeta, sigma, delta) + beta_values

t = np.linspace(0, 365*3, 365*3)  # Time vector (3 years)


#%% ODE model
def model(y, t, params):
    # Unpack state variables and parameters
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y
    (birth_rate, death_rate, theta, p_mortality, p_recover, tau,
     phi_mortality, phi_recover, zeta, sigma, delta, beta_l, beta_h) = params

    # Equations for the model
    N = S + Eh + Indh + Idh + Rh + El + Indl + Idl + Rl # Total population

    # Force of infection
    B_h = beta_h * (Indh + Idh)
    B_l = beta_l * Indl

    # Differential equations
    dSdt = birth_rate * N - (B_h + B_l) * S - death_rate * S + delta * (Rh + Rl)
    dEhdt = (B_h * S) - tau *(theta + (1-theta))*Eh
    dEldt = (B_l * S) - tau *(theta + (1-theta))*El
    dIdhdt = tau * theta * Eh - death_rate * Idh - phi_mortality * p_mortality * zeta * Idh - phi_recover * p_recover * sigma * Idh
    dIdldt = tau *theta * El - death_rate * Idl - p_mortality * zeta * Idl - p_recover * sigma * Idl
    dIndhdt = tau *(1-theta) * Eh - death_rate * Indh - phi_mortality * zeta * Indh - phi_recover * sigma * Indh
    dIndldt = tau *(1-theta) * El - death_rate * Indl  -zeta * Indl - sigma * Indl
    dRhdt = phi_recover * sigma * (p_recover * Idh + Indh) - (death_rate + delta) * Rh
    dRldt = sigma * (p_recover * Idl + Indl) - (death_rate + delta) * Rl

    return dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt

#%% Solve ODEs
solution = odeint(model, pop_values, t, args=(parameters,))

#%% Collect results in a DataFrame
columns = ['Susceptible', 'Exposed_High', 'Infected_NotDrug_High', 'Infected_Drug_High', 'Recovered_High',
           'Exposed_Low', 'Infected_NotDrug_Low', 'Infected_Drug_Low', 'Recovered_Low']
results = pd.DataFrame(solution, columns=columns)

Sdt, Ehdt, Indhdt, Idhdt, Rhdt, Eldt, Indldt, Idldt, Rldt = solution.T

fig = plt.figure(figsize=(12, 8), facecolor='white')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, Sdt, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, Ehdt, 'y', alpha=0.5, lw=2, label='Exposed High')
ax.plot(t, Indhdt, 'r', alpha=0.5, lw=2, label='Infected Not Drug High')
ax.plot(t, Idhdt, 'm', alpha=0.5, lw=2, label='Infected Drug High')
ax.plot(t, Rhdt, 'g', alpha=0.5, lw=2, label='Recovered High')
ax.plot(t, Eldt, 'c', alpha=0.5, lw=2, label='Exposed Low')
ax.plot(t, Indldt, color='orange', alpha=0.5, lw=2, label='Infected Not Drug Low')
ax.plot(t, Idldt, color='brown', alpha=0.5, lw=2, label='Infected Drug Low')
ax.plot(t, Rldt, color='olive', alpha=0.5, lw=2, label='Recovered Low')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number')
ax.set_ylim(0, 10000)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()

# %%
