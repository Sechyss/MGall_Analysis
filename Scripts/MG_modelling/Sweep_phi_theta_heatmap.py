# sweep_phi_theta_heatmap.py
# This script performs a 2D parameter sweep of phi_transmission vs theta for the SEIRS model,
# collects metrics (peak and equilibrium prevalence of the high-virulence strain) and
# produces heatmaps and a CSV table of results.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Models.SEIRS_Models import SEIRS_first_model
import pandas as pd

#%% --- Parameter Sweep: phi_transmission vs theta ---
# --- baseline parameters ---
# beta_l: baseline transmission rate for the low-virulence strain
beta_l = 0.25
# birth_rate / death_rate: demographic turnover; set to zero for short-term epidemic dynamics
birth_rate = 0.0
death_rate = 0.0
# delta: rate of waning immunity (1 / duration of immunity in days)
delta = 1/90
# delta_d: rate of starting drug (1 / average delay to drug start in days)
delta_d = 1/3
# p_recover: effectiveness of the drug in promoting recovery (fraction/scale)
p_recover = 0.5
# phi_recover: modifier for recovery in the high-virulence strain (how virulence affects recovery)
phi_recover = 0.75
# sigma: recovery rate (1 / infectious period)
sigma = 1/10
# tau: progression rate from exposed to infectious (1 / incubation period)
tau = 1/3

# --- initial conditions (normalized) ---
# Start with large susceptible population and a few initial infected for both strains.
S0, Eh0, Indh0, Idh0, Rh0, El0, Indl0, Idl0, Rl0 = 10000, 0, 5, 0, 0, 0, 5, 0, 0
y0 = np.array([S0, Eh0, Indh0, Idh0, Rh0, El0, Indl0, Idl0, Rl0])
# Normalize to proportions (model expects fractions/proportions rather than absolute counts)
y0 = y0 / y0.sum()

# --- time vector ---
# Simulate for one year with daily resolution
t = np.linspace(0, 365, 365)

# --- parameter grids ---
# phi_vals: grid for phi_transmission, the multiplicative factor that increases transmission of high-virulence strain
phi_vals = np.linspace(0.5, 3.0, 25)    # horizontal axis in heatmaps
# theta_vals: grid for theta, the fraction of exposed individuals who will take the drug
theta_vals = np.linspace(0, 0.9, 20)    # vertical axis in heatmaps

# storage arrays: shape (len(theta_vals), len(phi_vals))
# peak_high: maximum prevalence of high-virulence infection over the simulation (peak)
peak_high = np.zeros((len(theta_vals), len(phi_vals)))
# eq_high_abs: mean absolute prevalence of high-virulence infection near the end of the simulation (approx equilibrium)
eq_high_abs = np.zeros_like(peak_high)
# eq_high_frac: fraction of total infections that are high-virulence at equilibrium
eq_high_frac = np.zeros_like(peak_high)

# --- main parameter sweep ---
# Loop over theta (rows) and phi (columns). For every pair, run the ODE model and compute metrics.
for i, theta in enumerate(theta_vals):
    for j, phi in enumerate(phi_vals):
        # Pack parameters to match SEIRS_first_model signature:
        # (beta_l, birth_rate, death_rate, delta, delta_d, p_recover, phi_recover, phi_transmission, sigma, tau, theta)
        params = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
                  phi_recover, phi, sigma, tau, theta)

        # Integrate the ODE system for the given parameter set
        sol = odeint(SEIRS_first_model, y0, t, args=(params,))

        # Unpack solution columns into named variables for clarity
        S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = sol.T

        # total_high: combined prevalence of high-virulence infecteds (both drug and not-on-drug compartments)
        total_high = Indh + Idh
        # total_low: combined prevalence of low-virulence infecteds
        total_low = Indl + Idl
        # total_inf: total infected prevalence (high + low)
        total_inf = total_high + total_low

        # Compute metrics:
        # 1) peak_high: maximum value of total_high during the simulation
        peak_high[i, j] = total_high.max()

        # 2) eq_high_abs: approximate equilibrium absolute prevalence of high strain
        #    Use last 30 timepoints (last ~30 days) mean as a crude equilibrium estimate
        eq_high_abs[i, j] = total_high[-30:].mean()

        # 3) eq_high_frac: fraction of infections due to high strain at equilibrium
        #    add small epsilon in denominator to avoid division by zero
        eq_high_frac[i, j] = eq_high_abs[i, j] / (total_inf[-30:].mean() + 1e-12)

# --- plotting ---
# Create 3 side-by-side heatmaps: peak, equilibrium absolute, equilibrium fraction
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

def plot_heatmap(ax, data, title, cmap='viridis'):
    """
    Helper to draw a heatmap with proper axis labels and colorbar.

    - ax: matplotlib Axes
    - data: 2D numpy array with shape (len(theta_vals), len(phi_vals))
    - title: plot title
    - cmap: colormap
    """
    # imshow with origin='lower' so that smaller theta is at the bottom (natural orientation)
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=cmap,
                   extent=[phi_vals[0], phi_vals[-1], theta_vals[0], theta_vals[-1]])
    # Axis labels describe what each axis means
    ax.set_xlabel('phi_transmission (virulence transmissibility multiplier)')
    ax.set_ylabel('theta (fraction using symptom-blocking drug)')
    ax.set_title(title)
    # Add a small colorbar to the axis
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# Plot each metric with a visually distinct colormap
plot_heatmap(axes[0], peak_high, 'Peak high-virulence infection', cmap='plasma')
plot_heatmap(axes[1], eq_high_abs, 'Equilibrium high-virulence infection', cmap='magma')
plot_heatmap(axes[2], eq_high_frac, 'High strain fraction at equilibrium', cmap='inferno')

plt.tight_layout()
# Save figure to disk (Figures folder expected to exist)
plt.savefig('./Figures/phi_theta_heatmap_firstmodel.png', dpi=300)
# Show interactively if possible (depending on backend); for scripts this may save only
plt.show()

# --- optional: export results as table ---
# Flatten grids into long-form table for downstream analysis / plotting in R / Excel
df = pd.DataFrame({
    'theta': np.repeat(theta_vals, len(phi_vals)),               # repeat theta for each phi
    'phi_transmission': np.tile(phi_vals, len(theta_vals)),      # tile phi across theta rows
    'peak_high': peak_high.flatten(),                            # flatten 2D -> 1D
    'eq_high_abs': eq_high_abs.flatten(),
    'eq_high_frac': eq_high_frac.flatten()
})
# Save numeric results to a CSV table for reproducibility
df.to_csv('./Tables/phi_theta_heatmap_results_firstmodel.csv', index=False)
print("Saved numeric results to: ./Tables/phi_theta_heatmap_results_firstmodel.csv")
