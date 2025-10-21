# sweep_phi_theta_heatmap.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

# --- your model function ---
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
        sol = odeint(model, y0, t, args=(params,))
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
plt.savefig('phi_theta_heatmap.png', dpi=300)
plt.show()

# --- optional: export results ---
df = pd.DataFrame({
    'theta': np.repeat(theta_vals, len(phi_vals)),
    'phi_transmission': np.tile(phi_vals, len(theta_vals)),
    'peak_high': peak_high.flatten(),
    'eq_high_abs': eq_high_abs.flatten(),
    'eq_high_frac': eq_high_frac.flatten()
})
df.to_csv('phi_theta_heatmap_results.csv', index=False)
print("Saved numeric results to: phi_theta_heatmap_results.csv")
