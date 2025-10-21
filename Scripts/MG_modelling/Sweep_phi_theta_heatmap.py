# sweep_phi_theta_heatmap.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Models.SEIRS_Models import SEIRS_first_model
import pandas as pd
#%% --- Parameter Sweep: phi_transmission vs theta ---
# --- baseline parameters ---
beta_l = 0.25
birth_rate = 0.0
death_rate = 0.0
delta = 1/90
delta_d = 1/3
p_recover = 0.5
phi_recover = 0.75
sigma = 1/10
tau = 1/3

# --- initial conditions (normalized) ---
S0, Eh0, Indh0, Idh0, Rh0, El0, Indl0, Idl0, Rl0 = 10000, 0, 5, 0, 0, 0, 5, 0, 0
y0 = np.array([S0, Eh0, Indh0, Idh0, Rh0, El0, Indl0, Idl0, Rl0])
y0 = y0 / y0.sum()

# --- time vector ---
t = np.linspace(0, 365, 365)

# --- parameter grids ---
phi_vals = np.linspace(0.5, 3.0, 25)    # horizontal axis
theta_vals = np.linspace(0, 0.9, 20)    # vertical axis

# storage arrays
peak_high = np.zeros((len(theta_vals), len(phi_vals)))
eq_high_frac = np.zeros_like(peak_high)
eq_high_abs = np.zeros_like(peak_high)

# --- main parameter sweep ---
for i, theta in enumerate(theta_vals):
    for j, phi in enumerate(phi_vals):
        params = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
                  phi_recover, phi, sigma, tau, theta)
        sol = odeint(SEIRS_first_model, y0, t, args=(params,))
        S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = sol.T

        total_high = Indh + Idh
        total_low = Indl + Idl
        total_inf = total_high + total_low

        peak_high[i, j] = total_high.max()
        eq_high_abs[i, j] = total_high[-30:].mean()
        eq_high_frac[i, j] = eq_high_abs[i, j] / (total_inf[-30:].mean() + 1e-12)

# --- plotting ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

def plot_heatmap(ax, data, title, cmap='viridis'):
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=cmap,
                   extent=[phi_vals[0], phi_vals[-1], theta_vals[0], theta_vals[-1]])
    ax.set_xlabel('phi_transmission (virulence transmissibility multiplier)')
    ax.set_ylabel('theta (fraction using symptom-blocking drug)')
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

plot_heatmap(axes[0], peak_high, 'Peak high-virulence infection', cmap='plasma')
plot_heatmap(axes[1], eq_high_abs, 'Equilibrium high-virulence infection', cmap='magma')
plot_heatmap(axes[2], eq_high_frac, 'High strain fraction at equilibrium', cmap='inferno')

plt.tight_layout()
plt.savefig('./Figures/phi_theta_heatmap.png', dpi=300)
plt.show()

# --- optional: export results ---
df = pd.DataFrame({
    'theta': np.repeat(theta_vals, len(phi_vals)),
    'phi_transmission': np.tile(phi_vals, len(theta_vals)),
    'peak_high': peak_high.flatten(),
    'eq_high_abs': eq_high_abs.flatten(),
    'eq_high_frac': eq_high_frac.flatten()
})
df.to_csv('./Tables/phi_theta_heatmap_results.csv', index=False)
print("Saved numeric results to: ./Tables/phi_theta_heatmap_results.csv")
