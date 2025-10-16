#%% Imports
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#%% Initial Conditions (Proportions)
S = 10000
Eh = 0
Indh = 5
Idh = 0
Rh = 0
El = 0
Indl = 5
Idl = 0
Rl = 0

pop_values = np.array([S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl])
pop_values = pop_values / np.sum(pop_values)  # Normalize to proportions

#%% Parameters
theta = 0.6  # Proportion of exposed who will eventually take the drug
p_recover = 1.0  # Drug effect on recovery
phi_transmission = 2.0  # High virulence transmission multiplier
phi_recover = 0.75  # High virulence recovery reduction
sigma = 1/7  # Recovery rate
delta = 1/90  # Immunity loss rate
tau = 1/5  # Progression from exposed to infectious
delta_d = 1/3  # Delay rate for starting drug (~3 days)

birth_rate = 0.0
death_rate = 0.0

beta_l = 0.5
beta_h = phi_transmission * beta_l

parameters = (birth_rate, death_rate, theta, p_recover, tau,
              phi_recover, sigma, delta, beta_l, beta_h, delta_d)

t = np.linspace(0, 365*1, 365*1)  # 1 year, daily resolution

#%% ODE Model
def model(y, t, params):
    # Prevent negative values
    y = [max(0, val) for val in y]
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y

    # Unpack parameters
    (birth_rate, death_rate, theta, p_recover, tau,
     phi_recover, sigma, delta, beta_l, beta_h, delta_d) = params

    # Force of infection (proportions)
    B_h = beta_h * (Indh + Idh)
    B_l = beta_l * (Indl + Idl)

    # Differential equations
    dSdt = birth_rate - (B_h + B_l) * S + delta * (Rh + Rl)
    dEhdt = B_h * S - tau * Eh
    dEldt = B_l * S - tau * El

    # High virulence infected dynamics
    dIndhdt = tau * Eh - delta_d * theta * Indh - phi_recover * sigma * Indh
    dIdhdt = delta_d * theta * Indh - phi_recover * p_recover * sigma * Idh

    # Low virulence infected dynamics
    dIndldt = tau * El - delta_d * theta * Indl - sigma * Indl
    dIdldt = delta_d * theta * Indl - p_recover * sigma * Idl

    # Recovery dynamics
    dRhdt = phi_recover * sigma * (p_recover * Idh + Indh * (1 - theta)) - delta * Rh
    dRldt = sigma * (p_recover * Idl + Indl * (1 - theta)) - delta * Rl

    return dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt

#%% Solve ODEs
solution = odeint(model, pop_values, t, args=(parameters,))

#%% Collect results
columns = ['Susceptible', 'Exposed_High', 'Infected_NotDrug_High', 'Infected_Drug_High', 'Recovered_High',
           'Exposed_Low', 'Infected_NotDrug_Low', 'Infected_Drug_Low', 'Recovered_Low']
results = pd.DataFrame(solution, columns=columns)

Sdt, Ehdt, Indhdt, Idhdt, Rhdt, Eldt, Indldt, Idldt, Rldt = solution.T

#%% Plotting
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
ax.set_ylabel('Proportion of Population')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.tight_layout()
plt.show()