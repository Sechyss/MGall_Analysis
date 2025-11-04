import numpy as np

#%% ODE Model
def SEIRS_first_model(y, t, params):
    """
    SEIRS_first_model
    - State vector y (length 9): [S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl]
    - params (tuple length 11): (beta_l, birth_rate, death_rate, delta, delta_d,
                                p_recover, phi_recover, phi_transmission, sigma, tau, theta)

    NOTE: This variant INTENTIONALLY ignores Idl when computing low-strain
    force-of-infection (B_l = beta_l * Indl). Use SEIRS_second_model if treated
    low cases (Idl) should contribute to transmission.

    Improvements:
    - use numpy arrays for numeric safety and performance
    - clip negative states to zero using np.maximum
    - validate params length and basic mass-balance sanity
    """
    # ensure numeric ndarray and non-negative compartments
    y = np.maximum(np.asarray(y, dtype=float), 0.0)
    S, Eh, Indh, Idh, Rh, El, Indl, Idl, Rl = y

    # basic validation to catch wrong parameter ordering / length early
    if not hasattr(params, "__len__") or len(params) != 11:
        raise ValueError("params must be a sequence of length 11: "
                         "(beta_l, birth_rate, death_rate, delta, delta_d, "
                         "p_recover, phi_recover, phi_transmission, sigma, tau, theta)")

    (beta_l, birth_rate, death_rate, delta, delta_d, p_recover,
     phi_recover, phi_transmission, sigma, tau, theta) = params

    # compute transmission contributions
    beta_h = phi_transmission * beta_l
    # high-strain: include treated and untreated high infecteds
    B_h = beta_h * (Indh + Idh)
    # low-strain: INTENTIONALLY only untreated low infecteds contribute in this variant
    B_l = beta_l * Indl

    # differential equations
    dSdt = birth_rate - (B_h + B_l) * S + delta * (Rh + Rl) - death_rate * S
    dEhdt = B_h * S - tau * Eh - death_rate * Eh
    dEldt = B_l * S - tau * El - death_rate * El

    dIndhdt = tau * Eh - delta_d * theta * Indh - phi_recover * sigma * Indh - death_rate * Indh
    dIdhdt = delta_d * theta * Indh - phi_recover * p_recover * sigma * Idh - death_rate * Idh

    dIndldt = tau * El - delta_d * theta * Indl - sigma * Indl - death_rate * Indl
    dIdldt = delta_d * theta * Indl - p_recover * sigma * Idl - death_rate * Idl

    dRhdt = phi_recover * sigma * (p_recover * Idh + Indh) - delta * Rh - death_rate * Rh
    dRldt = sigma * (p_recover * Idl + Indl) - delta * Rl - death_rate * Rl

    # quick mass-balance check (optional small tolerance)
    total = S + Eh + Indh + Idh + Rh + El + Indl + Idl + Rl
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("Non-finite or non-positive total population in SEIRS_first_model")
    # return as tuple (compatible with odeint)
    return dSdt, dEhdt, dIndhdt, dIdhdt, dRhdt, dEldt, dIndldt, dIdldt, dRldt
