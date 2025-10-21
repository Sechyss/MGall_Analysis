#%% Imports
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy.linalg import eigvals
from Models.SEIRS_Models import SEIRS_first_model

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
theta = 0.0  # Proportion taking the drug
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

#%% Solve ODEs
solution = odeint(SEIRS_first_model, pop_values, t, args=(parameters,))
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
plt.savefig('./Figures/model_deeper_dive_dynamics.png', dpi=300)
plt.show()

#%% 1️⃣ Check equilibrium
steady_state = results.iloc[-1]
print("\n--- Approximate Equilibrium (final day) ---")
print(steady_state)

#%% 2️⃣ Compute Basic Reproduction Numbers (approximation)
R0_low = beta_l / sigma
R0_high = (phi_transmission * beta_l) / (phi_recover * sigma)
print(f"\nEstimated R0_low = {R0_low:.2f}")
print(f"Estimated R0_high = {R0_high:.2f}")

#%% 3️⃣ Jacobian-based stability at DFE
def jacobian_at_dfe(params):
    (beta_l, br, dr, delta, delta_d, p_rec,
     phi_rec, phi_trans, sigma, tau, theta) = params
    beta_h = phi_trans * beta_l
    J = np.zeros((4, 4))
    # Rows/cols = [Eh, El, Indh, Indl]
    J[0, 0] = -tau - dr
    J[0, 2] = beta_h
    J[1, 1] = -tau - dr
    J[1, 3] = beta_l
    J[2, 0] = tau
    J[2, 2] = -(phi_rec * sigma + dr)
    J[3, 1] = tau
    J[3, 3] = -(sigma + dr)
    return J

J = jacobian_at_dfe(parameters)
eigs = eigvals(J)
print("\nJacobian eigenvalues at DFE:", np.round(eigs, 4))
print("Stable equilibrium?", np.all(np.real(eigs) < 0))

#%% 4️⃣ Phase plane: high vs low infection
plt.figure(figsize=(7, 6))
plt.plot(Indldt + Idldt, Indhdt + Idhdt, color='purple', lw=2)
plt.xlabel("Low-virulence infection proportion")
plt.ylabel("High-virulence infection proportion")
plt.title("Phase plane: competition between strains")
plt.grid(True, alpha=0.4)
plt.savefig('./Figures/phase_plane_high_low_infection.png', dpi=300)
plt.show()

#%% 5️⃣ Scenario comparison: vary theta (drug usage)
thetas = np.linspace(0, 0.9, 10)
peak_high, peak_low = [], []
eq_high, eq_low = [], []

for th in thetas:
    params_var = (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
                  phi_recover, phi_transmission, sigma, tau, th)
    sol = odeint(model, pop_values, t, args=(params_var,))
    _, _, Indh, Idh, _, _, Indl, Idl, _ = sol.T
    total_high = Indh + Idh
    total_low = Indl + Idl
    peak_high.append(np.max(total_high))
    peak_low.append(np.max(total_low))
    eq_high.append(np.mean(total_high[-30:]))  # equilibrium approx (last 30 days)
    eq_low.append(np.mean(total_low[-30:]))

#%% 6️⃣ Plot: effect of theta on infection peaks
plt.figure(figsize=(7, 5))
plt.plot(thetas, peak_low, 'o-', label='Low virulence peak')
plt.plot(thetas, peak_high, 's-', label='High virulence peak')
plt.xlabel("θ (fraction using symptom-blocking drug)")
plt.ylabel("Peak infection proportion")
plt.title("Effect of drug coverage on infection peaks")
plt.legend()
plt.tight_layout()
plt.savefig('./Figures/drug_coverage_infection_peaks.png', dpi=300)
plt.show()

#%% 7️⃣ Bifurcation-like plot: dominance at equilibrium
frac_high = np.array(eq_high) / (np.array(eq_high) + np.array(eq_low))
plt.figure(figsize=(7, 5))
plt.plot(thetas, frac_high, 'd-', color='crimson', lw=2)
plt.axhline(0.5, color='gray', ls='--', lw=1)
plt.xlabel("θ (fraction using symptom-blocking drug)")
plt.ylabel("High virulence fraction at equilibrium")
plt.title("Bifurcation-like diagram: selection for virulence vs drug coverage")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./Figures/bifurcation_virulence_drug_coverage.png', dpi=300)
plt.show()

print("\n--- Summary of equilibrium virulence dominance ---")
for th, f in zip(thetas, frac_high):
    dom = "High" if f > 0.5 else "Low"
    print(f"θ={th:.2f} → {f*100:.1f}% high-virulence (dominant: {dom})")
