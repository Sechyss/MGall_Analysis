#%% Imports
import numpy as np                # numerical arrays and math
import pandas as pd               # dataframes for results storage
from scipy.integrate import odeint # ODE integrator
import matplotlib.pyplot as plt   # plotting
import seaborn as sns             # nicer plotting styles
from Models.SEIRS_Models import SEIRS_first_model  # the ODE model function to integrate
from Models import params as model_params  # import default parameters

# Make any stochastic behaviour reproducible
np.random.seed(42)

#%% --- Load Initial Conditions from module ---
# --- baseline parameters (loaded from Models.params, with sensible defaults) ---
# use getattr so script still runs if a name is missing in model_params
beta_l = getattr(model_params, "beta_l", 0.25)
birth_rate = getattr(model_params, "birth_rate", 0.0)
death_rate = getattr(model_params, "death_rate", 0.0)
delta = getattr(model_params, "delta", 1/90)
delta_d = getattr(model_params, "delta_d", 1/3)
p_recover = getattr(model_params, "p_recover", 0.5)
phi_recover = getattr(model_params, "phi_recover", 0.75)
sigma = getattr(model_params, "sigma", 1/10)
tau = getattr(model_params, "tau", 1/3)
# optional defaults for grids / time (override in params.py if desired)
phi_transmission = getattr(model_params, "phi_transmission", None)
theta = getattr(model_params, "theta", None)
t_max = getattr(model_params, "t_max", 365)
t_steps = getattr(model_params, "t_steps", 365)

# --- initial conditions (loaded from Models.params) ---
# load initial compartment counts from params.py (fall back to previous defaults if missing)
S0 = getattr(model_params, "S", 10000)
Eh0 = getattr(model_params, "Eh", 0)
Indh0 = getattr(model_params, "Indh", 5)
Idh0 = getattr(model_params, "Idh", 0)
Rh0 = getattr(model_params, "Rh", 0)
El0 = getattr(model_params, "El", 0)
Indl0 = getattr(model_params, "Indl", 5)
Idl0 = getattr(model_params, "Idl", 0)
Rl0 = getattr(model_params, "Rl", 0)
y0 = np.array([S0, Eh0, Indh0, Idh0, Rh0, El0, Indl0, Idl0, Rl0])
# Normalize to proportions (model expects fractions/proportions rather than absolute counts)
y0 = y0 / y0.sum()

parameters = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
              phi_recover, phi_transmission, sigma, tau, theta)
# --- time vector ---
# Use t_max / t_steps from params.py when available, otherwise daily for one year
t = np.linspace(0, t_max, int(t_steps))

#%% Solve ODEs (baseline run)
# Integrate the ODE system once at baseline parameters to produce a results DataFrame
solution = odeint(SEIRS_first_model, y0, t, args=(parameters,))

# Column names corresponding to model state ordering returned by SEIRS_first_model
columns = [
    'Susceptible',
    'Exposed_High',
    'Infected_NotDrug_High',
    'Infected_Drug_High',
    'Recovered_High',
    'Exposed_Low',
    'Infected_NotDrug_Low',
    'Infected_Drug_Low',
    'Recovered_Low'
]

# Store solution in a DataFrame for easy downstream analysis and plotting
results = pd.DataFrame(solution, columns=columns)

# Unpack solution columns into named arrays for plotting convenience
Sdt, Ehdt, Indhdt, Idhdt, Rhdt, Eldt, Indldt, Idldt, Rldt = solution.T

#%% --- Sensitivity Analysis Section ---

def outcome(df):
    """
    Define the scalar model outcome used for sensitivity analysis.
    Here: peak total infection (max across time of all infected compartments).
    """
    total_infected = (
        df['Infected_NotDrug_High'] + df['Infected_Drug_High'] +
        df['Infected_NotDrug_Low'] + df['Infected_Drug_Low']
    )
    return total_infected.max()

def run_model(param_dict):
    """
    Run the ODE model with parameters provided in param_dict and return the scalar outcome.
    - Build the params tuple in the same order SEIRS_first_model expects.
    - Integrate ODE and compute the chosen outcome from the resulting DataFrame.
    """
    # Unpack parameters from input dictionary into a tuple in the correct order
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

    # IMPORTANT: use the imported model function (SEIRS_first_model).
    # Original code used 'model' which is undefined here — that is a bug.
    sol = odeint(SEIRS_first_model, y0, t, args=(params,))

    # Convert to DataFrame using the same columns defined earlier
    df = pd.DataFrame(sol, columns=columns)

    # Return scalar metric for sensitivity analysis
    return outcome(df)

# Baseline parameters as a dictionary used for sensitivity reference
baseline = {
    'beta_l': 0.25,
    'birth_rate': 0.0,
    'death_rate': 0.0,
    'delta': 1/90,
    'delta_d': 1/3,
    'p_recover': 0.5,
    'phi_recover': 0.75,
    'phi_transmission': 1.3,
    'sigma': 1/10,
    'tau': 1/3,
    'theta': 0.3
}

def sensitivity_analysis(param_dict, variation=0.1):
    """
    One-at-a-time (OAT) sensitivity analysis.
    For each parameter in param_dict:
      - evaluate outcome at (1 - variation) * param and (1 + variation) * param
      - compute a simple sensitivity index = (high - low) / (2 * variation * baseline_value)
    Returns a DataFrame with low/high/baseline and sensitivity_index for each parameter.
    """
    baseline_result = run_model(param_dict)  # run model at baseline for reference
    sens_results = {}

    # Iterate through each parameter to perturb it independently
    for key in param_dict:
        # low and high perturbed values (±variation fraction)
        p_values = [param_dict[key] * (1 - variation), param_dict[key] * (1 + variation)]
        outcomes = []
        for val in p_values:
            test_params = param_dict.copy()
            test_params[key] = val
            outcomes.append(run_model(test_params))  # run model and store outcome

        # Calculate a normalized sensitivity index (finite-difference / fractional change)
        sens_results[key] = {
            'low': outcomes[0],
            'high': outcomes[1],
            'baseline': baseline_result,
            # avoid division by zero; index approximates derivative normalized by baseline
            'sensitivity_index': (outcomes[1] - outcomes[0]) / (2 * variation * param_dict[key])
        }

    # Return results as a DataFrame with parameters as rows
    return pd.DataFrame(sens_results).T

# Run the sensitivity analysis with default 10% perturbation
sens_df = sensitivity_analysis(baseline, variation=0.1)

# Print concise summary table showing baseline, low, high and sensitivity index
print("\n=== Sensitivity Analysis Results ===")
print(sens_df[['baseline', 'low', 'high', 'sensitivity_index']])

#%% --- Plot sensitivity results as a horizontal bar chart (beautified) ---

# Create figure and set seaborn theme for nicer style
plt.figure(figsize=(8, 6))
sns.set_theme(style="whitegrid")

# Sort parameters by absolute sensitivity (descending) to make the tornado plot readable
sens_sorted = sens_df.sort_values(by='sensitivity_index', key=np.abs, ascending=False)

# Draw horizontal bars centered at zero; use a diverging colormap mapped to number of params
plt.barh(
    sens_sorted.index,
    sens_sorted['sensitivity_index'],
    color=sns.color_palette('coolwarm', len(sens_sorted)),  # color for each bar
    edgecolor='black'
)

# Add vertical zero line for reference (positive effect vs negative effect)
plt.axvline(0, color='black', lw=1)
plt.xlabel('Sensitivity Index (Effect on Peak Infection)', fontsize=13)
plt.ylabel('Parameter', fontsize=13)
plt.title('Tornado Plot of Parameter Sensitivity', fontsize=15, weight='bold')

# Annotate each bar with its numeric sensitivity value
for i, (idx, val) in enumerate(zip(sens_sorted.index, sens_sorted['sensitivity_index'])):
    # place label slightly offset from bar end; choose alignment based on sign
    plt.text(val + np.sign(val) * 0.01, i, f"{val:.2f}", va='center',
             ha='left' if val >= 0 else 'right', fontsize=10)

plt.tight_layout()

# Save plot to file for reproducibility; plt.show() will also display if backend is interactive
plt.savefig('sensitivity_tornado_plot.png', dpi=300)
plt.show()

# %%
# End of script
