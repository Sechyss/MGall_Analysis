# Initial absolute counts for each compartment (counts)
S = 10000
Eh = 0
Indh = 5
Idh = 0
Rh = 0
El = 0
Indl = 5
Idl = 0
Rl = 0

# Model parameters (units: days for rates unless noted)
# Treatment / biology
theta = 0.0            # fraction of exposed who will eventually take the drug (0..1)
p_recover = 0.5        # effectiveness factor of drug on recovery (0..1)

# Virulence modifiers
phi_transmission = 1.3 # multiplier for high-virulence transmission (beta_h = phi_transmission * beta_l)
phi_recover = 0.75     # modifier for recovery for high-virulence (if <1, recovery is slower)

# Rates (per day)
sigma = 1/10           # recovery rate (1 / infectious period)
delta = 1/90           # immunity waning rate (1 / duration of immunity)
tau = 1/3              # progression from exposed -> infectious (1 / incubation period)
delta_d = 1/3          # rate of starting drug (1 / mean delay to drug start)

# Demography and transmission
birth_rate = 0.0
death_rate = 0.0
beta_l = 0.25          # baseline transmission rate for low-virulence strain

# Time grid defaults (used if scripts import these)
t_max = 365            # simulation length in days
t_steps = 365          # number of time points (e.g. daily resolution)

# Optional: recommended ranges (for reference)
# theta: 0.0 - 0.9
# p_recover: 0.1 - 0.95
# delta: 1/365 (1 year) - 1/30 (1 month)
