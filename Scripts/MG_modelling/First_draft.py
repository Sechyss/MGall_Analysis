#%% Imports
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy.linalg import eigvals

np.random.seed(42)

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
theta = 0.0  # Proportion of exposed who will eventually take the drug
p_recover = 0.5  # Drug effect on recovery
phi_transmission = 1.3  # High virulence transmission multiplier
phi_recover = 0.75  # High virulence recovery reduction
sigma = 1/10  # Recovery rate
delta = 1/90  # Immunity loss rate
tau = 1/3  # Progression from exposed to infectious
delta_d = 1/3  # Delay rate for starting drug (~3 days)
birth_rate = 0.0
death_rate = 0.0
beta_l = 0.25

parameters = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
              phi_recover, phi_transmission, sigma, tau, theta)

# Time vector: 1 year, daily resolution
t = np.linspace(0, 365, 365)

#%% ODE Model
def model(y, t, params):
    y = [max(0, val) for val in y]
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y

    (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
     phi_recover, phi_transmission, sigma, tau, theta) = params

    beta_h = phi_transmission * beta_l
    B_h = beta_h * (Indh + Idh)
    B_l = beta_l * Indl

    dSdt = birth_rate - (B_h + B_l) * S + delta * (Rh + Rl) - death_rate * S
    dEhdt = B_h * S - tau * Eh - death_rate * Eh
    dEldt = B_l * S - tau * El - death_rate * El

    dIndhdt = tau * Eh - delta_d * theta * Indh - phi_recover * sigma * Indh - death_rate * Indh
    dIdhdt = delta_d * theta * Indh - phi_recover * p_recover * sigma * Idh - death_rate * Idh

    dIndldt = tau * El - delta_d * theta * Indl - sigma * Indl - death_rate * Indl
    dIdldt = delta_d * theta * Indl - p_recover * sigma * Idl - death_rate * Idl

    dRhdt = phi_recover * sigma * (p_recover * Idh + Indh) - delta * Rh - death_rate * Rh
    dRldt = sigma * (p_recover * Idl + Indl) - delta * Rl - death_rate * Rl

    return dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt

#%% Solve ODEs
solution = odeint(model, pop_values, t, args=(parameters,))
columns = ['Susceptible', 'Exposed_High', 'Infected_NotDrug_High', 'Infected_Drug_High', 'Recovered_High',
           'Exposed_Low', 'Infected_NotDrug_Low', 'Infected_Drug_Low', 'Recovered_Low']
results = pd.DataFrame(solution, columns=columns)
Sdt, Ehdt, Indhdt, Idhdt, Rhdt, Eldt, Indldt, Idldt, Rldt = solution.T

#%% Plot time dynamics
fig = plt.figure(figsize=(12, 8), facecolor='white')
ax = fig.add_subplot(111, facecolor='#f4f4f4', axisbelow=True)
ax.plot(t, Sdt, 'b', lw=2, label='Susceptible')
ax.plot(t, Ehdt, 'y', lw=2, label='Exposed High')
ax.plot(t, Indhdt, 'r', lw=2, label='Infected Not Drug High')
ax.plot(t, Idhdt, 'm', lw=2, label='Infected Drug High')
ax.plot(t, Rhdt, 'g', lw=2, label='Recovered High')
ax.plot(t, Eldt, 'c', lw=2, label='Exposed Low')
ax.plot(t, Indldt, color='orange', lw=2, label='Infected Not Drug Low')
ax.plot(t, Idldt, color='brown', lw=2, label='Infected Drug Low')
ax.plot(t, Rldt, color='olive', lw=2, label='Recovered Low')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Proportion of Population')
ax.legend(framealpha=0.7)
plt.tight_layout()
plt.savefig('first_draft_model_dynamics.png', dpi=300)