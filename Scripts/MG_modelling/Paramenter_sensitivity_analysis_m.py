#%% Imports
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
from Models import SEIRS_first_model

np.random.seed(42)

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
theta = 0.3  # Proportion of exposed who will eventually take the drug
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


# Pack parameters into a tuple
parameters = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
     phi_recover, phi_transmission, sigma, tau, theta)

# Time vector: 1 year, daily resolution
t = np.linspace(0, 365*1, 365*1)

#%% Solve ODEs
# Integrate the ODE system
solution = odeint(SEIRS_first_model, pop_values, t, args=(parameters,))

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
        param_dict['beta_l'],
        param_dict['birth_rate'],
        param_dict['death_rate'],
        param_dict['delta'],
        param_dict['delta_d'],
        param_dict['p_recover'],
        param_dict['phi_recover'],
        param_dict['phi_transmission'],
        param_dict['sigma'],
        param_dict['tau'],
        param_dict['theta']
    )
    sol = odeint(model, pop_values, t, args=(params,))
    df = pd.DataFrame(sol, columns=columns)
    return outcome(df)

# Baseline parameters as dictionary
baseline = {
    'beta_l': 0.25,
    'birth_rate': 0.0,
    'death_rate': 0.0,
    'delta': 1/90,
    'delta_d': 1/3,
    'p_recover': 0.5,
    'phi_recover': 0.75,
    'phi_transmission': 1.5,
    'sigma': 1/10,
    'tau': 1/3,
    'theta': 0.3
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

plt.figure(figsize=(8, 6))
sns.set_theme(style="whitegrid")

# Sort parameters by absolute sensitivity (descending)
sens_sorted = sens_df.sort_values(by='sensitivity_index', key=np.abs, ascending=False)

# Plot horizontal bars centered at zero
plt.barh(
    sens_sorted.index,
    sens_sorted['sensitivity_index'],
    color=sns.color_palette('coolwarm', len(sens_sorted)),
    edgecolor='black'
)

plt.axvline(0, color='black', lw=1)
plt.xlabel('Sensitivity Index (Effect on Peak Infection)', fontsize=13)
plt.ylabel('Parameter', fontsize=13)
plt.title('Tornado Plot of Parameter Sensitivity', fontsize=15, weight='bold')

# Annotate values
for i, (idx, val) in enumerate(zip(sens_sorted.index, sens_sorted['sensitivity_index'])):
    plt.text(val + np.sign(val)*0.01, i, f"{val:.2f}", va='center',
             ha='left' if val >= 0 else 'right', fontsize=10)

plt.tight_layout()
plt.savefig('sensitivity_tornado_plot.png', dpi=300)
plt.show()

# %%
