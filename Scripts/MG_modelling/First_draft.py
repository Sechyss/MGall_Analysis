import numpy as np
import pandas as pd
from scipy.integrate import odeint

#%% Parameters

# Initial population values
S = 2000  # Initial susceptible population
Eh = 0  # Initial exposed humans high virulence pathogen
Indh = 0  # Initial infected humans high virulence pathogen not taking the drug
Idh = 0  # Initial infected humans high virulence pathogen taking the drug
Rh = 0  # Initial recovered humans high virulence pathogen
El = 0  # Initial exposed humans low virulence pathogen
Indl = 0  # Initial infected humans low virulence pathogen not taking the drug
Idl = 0  # Initial infected humans low virulence pathogen taking the drug
Rl = 0  # Initial recovered humans low virulence pathogen

pop_values = (S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl)
N = sum(pop_values)  # Total initial population


# Model parameters
birth_rate = 0.02 # Birth rate
death_rate = 0.01 # Natural death rate due to other causes

theta = 0.1  # Proportion of population that takes the drug
p_mortality = 0.1  # Protection effect of the drug against mortality
p_recover = 0.1  # Protection effect of the drug for recovery
lambda_transmission = 1.5  # Increase in transmission rate for high virulence pathogen
lambda_mortality = 2  # Increase in mortality rate for high virulence pathogen
lambda_recover = 0.5  # Decrease in recovery rate for high virulence pathogen
zeta = 0.1  # Death rate due to infection
sigma = 1/7  # Rate of recovery from infection
delta = 1/150  # Rate of loss of immunity

beta_l = 0.3  # Transmission rate low virulence pathogen
beta_h = lambda_transmission * beta_l  # Transmission rate high virulence pathogen

beta_values = beta_l, beta_h
parameters = (birth_rate, death_rate, theta, p_mortality, p_recover, lambda_transmission,
              lambda_mortality, lambda_recover, zeta, sigma, delta) + beta_values

t = np.linspace(0, 365*3, 365*3)  # Time vector (3 years)


#%% ODE model
def model(y, t, params):
    # Unpack state variables and parameters
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y
    (birth_rate, death_rate, theta, p_mortality, p_recover, lambda_transmission,
     lambda_mortality, lambda_recover, zeta, sigma, delta, beta_l, beta_h) = params

    # Equations for the model
    N = S + Eh + Indh + Idh + Rh + El + Indl + Idl + Rl # Total population

    # Force of infection
    B_h = beta_h * (Indh + Idh)
    B_l = beta_l * Indl

    # Differential equations
    dSdt = birth_rate * N - (B_h + B_l) * S - death_rate * S + delta * (Rh + Rl)
    dEhdt = (B_h * S) - (theta + (1-theta))*Eh
    dEldt = (B_l * S) - (theta + (1-theta))*El
    dIdhdt = theta * Eh - death_rate * Idh - lambda_mortality * p_mortality * zeta * Idh - lambda_recover * p_recover * sigma * Idh
    dIdldt = theta * El - death_rate * Idl - p_mortality * zeta * Idl - p_recover * sigma * Idl
    dIndhdt = (1-theta) * Eh - death_rate * Indh - lambda_mortality * zeta * Indh - lambda_recover * sigma * Indh
    dIndldt = (1-theta) * El - death_rate * Indl  -zeta * Indl - sigma * Indl
    dRhdt = lambda_recover * sigma * (p_recover * Idh + Indh) - (death_rate + delta) * Rh
    dRldt = sigma * (p_recover * Idl + Indl) - (death_rate + delta) * Rl

    return [dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt]

#%% Solve ODEs
solution = odeint(model, pop_values, t, args=(parameters,))

# Collect results in a DataFrame

columns = ['Susceptible', 'Exposed_High', 'Infected_NotDrug_High', 'Infected_Drug_High', 'Recovered_High',
           'Exposed_Low', 'Infected_NotDrug_Low', 'Infected_Drug_Low', 'Recovered_Low']
results = pd.DataFrame(solution, columns=columns)