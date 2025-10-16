#%% Imports
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns

#%% Initial Conditions (Proportions)
# Set initial population values for each compartment
S = 10000
Eh = 0
Indh = 5
Idh = 0
Rh = 0
El = 0
Indl = 5
Idl = 0
Rl = 0

# Normalize to proportions
pop_values = np.array([S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl])
pop_values = pop_values / np.sum(pop_values)

#%% Parameters
# Define model parameters
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

# Pack parameters into a tuple
parameters = (birth_rate, death_rate, theta, p_recover, tau,
              phi_recover, sigma, delta, beta_l, beta_h, delta_d)

# Time vector: 1 year, daily resolution
t = np.linspace(0, 365*1, 365*1)

#%% ODE Model
def model(y, t, params):
    # Prevent negative values in compartments
    y = [max(0, val) for val in y]
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y

    (birth_rate, death_rate, theta, p_recover, tau,
     phi_recover, sigma, delta, beta_l, beta_h, delta_d) = params

    # Calculate force of infection for high and low virulence
    B_h = beta_h * Indh
    B_l = beta_l * (Indl + Idl)

    # Differential equations for each compartment
    dSdt = birth_rate - (B_h + B_l) * S + delta * (Rh + Rl)
    dEhdt = B_h * S - tau * Eh
    dEldt = B_l * S - tau * El

    # High virulence infected
    dIndhdt = tau * Eh - delta_d * theta * Indh - phi_recover * sigma * Indh
    dIdhdt = delta_d * theta * Indh - phi_recover * p_recover * sigma * Idh

    # Low virulence infected
    dIndldt = tau * El - delta_d * theta * Indl - sigma * Indl
    dIdldt = delta_d * theta * Indl - p_recover * sigma * Idl

    # Recovery compartments
    dRhdt = phi_recover * sigma * (p_recover * Idh + Indh * (1 - theta)) - delta * Rh
    dRldt = sigma * (p_recover * Idl + Indl * (1 - theta)) - delta * Rl

    return dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt

#%% Solve ODEs
# Integrate the ODE system
solution = odeint(model, pop_values, t, args=(parameters,))

# Create DataFrame for results
columns = ['Susceptible', 'Exposed_High', 'Infected_NotDrug_High', 'Infected_Drug_High', 'Recovered_High',
           'Exposed_Low', 'Infected_NotDrug_Low', 'Infected_Drug_Low', 'Recovered_Low']
results = pd.DataFrame(solution, columns=columns)
Sdt, Ehdt, Indhdt, Idhdt, Rhdt, Eldt, Indldt, Idldt, Rldt = solution.T

#%% --- Sensitivity Analysis Section ---

def outcome(df):
    """Define the model outcome metric (peak total infection)."""
    # Sum all infected compartments and return the maximum value (peak)
    total_infected = (df['Infected_NotDrug_High'] + df['Infected_Drug_High'] +
                      df['Infected_NotDrug_Low'] + df['Infected_Drug_Low'])
    return total_infected.max()

def run_model(param_dict):
    """Run the ODE model with the given parameter dictionary."""
    # Unpack parameters from dictionary and run the model
    params = (
        param_dict['birth_rate'], param_dict['death_rate'], param_dict['theta'],
        param_dict['p_recover'], param_dict['tau'], param_dict['phi_recover'],
        param_dict['sigma'], param_dict['delta'], param_dict['beta_l'],
        param_dict['beta_h'], param_dict['delta_d']
    )
    sol = odeint(model, pop_values, t, args=(params,))
    df = pd.DataFrame(sol, columns=columns)
    return outcome(df)

# Baseline parameters as dictionary
baseline = {
    'birth_rate': 0.0, 'death_rate': 0.0, 'theta': 0.6, 'p_recover': 1.0, 'tau': 1/5,
    'phi_recover': 0.75, 'sigma': 1/7, 'delta': 1/90, 'beta_l': 0.5,
    'beta_h': 2.0 * 0.5, 'delta_d': 1/3
}

def sensitivity_analysis(param_dict, variation=0.1):
    """Run one-at-a-time sensitivity analysis.
    For each parameter, decrease and increase it by 'variation' percent and record the effect on the outcome."""
    baseline_result = run_model(param_dict)
    sens_results = {}
    for key in param_dict:
        # Test parameter at low and high values (±variation)
        p_values = [param_dict[key]*(1 - variation), param_dict[key]*(1 + variation)]
        outcomes = []
        for val in p_values:
            test_params = param_dict.copy()
            test_params[key] = val
            outcomes.append(run_model(test_params))
        # Store results and calculate a simple sensitivity index
        sens_results[key] = {
            'low': outcomes[0],
            'high': outcomes[1],
            'baseline': baseline_result,
            'sensitivity_index': (outcomes[1] - outcomes[0]) / (2 * variation * param_dict[key])
        }
    return pd.DataFrame(sens_results).T

# Run and display sensitivity results
sens_df = sensitivity_analysis(baseline, variation=0.1)
print("\n=== Sensitivity Analysis Results ===")
print(sens_df[['baseline', 'low', 'high', 'sensitivity_index']])

# Plot sensitivity results as a horizontal bar chart (beautified)


plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

# Sort by absolute sensitivity for better visual order
sens_df_sorted = sens_df.reindex(sens_df['sensitivity_index'].abs().sort_values(ascending=True).index)

barplot = sns.barplot(
    x='sensitivity_index',
    y=sens_df_sorted.index,
    data=sens_df_sorted,
    palette='viridis',
    orient='h'
)

# Add value labels to bars
for i, (value, param) in enumerate(zip(sens_df_sorted['sensitivity_index'], sens_df_sorted.index)):
    plt.text(
        value + 0.01 * np.sign(value),  # Offset label slightly from bar
        i,
        f"{value:.2f}",
        va='center',
        ha='left' if value >= 0 else 'right',
        fontsize=10
    )

plt.xlabel('Sensitivity Index', fontsize=14)
plt.ylabel('Parameter', fontsize=14)
plt.title('Parameter Sensitivity on Peak Infections', fontsize=16, weight='bold')
plt.tight_layout()
plt.show()
