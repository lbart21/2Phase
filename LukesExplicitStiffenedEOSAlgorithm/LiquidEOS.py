def update_from_rhoP(state, coefficients):
    p = state["p"]
    rho = state["rho"]
    gamma = coefficients["gamma"]
    Cv = coefficients["Cv"]
    Cp = coefficients["Cp"]
    p_inf = coefficients["p_inf"]
    b = coefficients["b"]
    q = coefficients["q"]
    T = (1 / rho - b) * (p + p_inf) / (Cv * (gamma - 1))
    state["T"] = T
    state["a"] = (gamma * 1 / rho ** 2 * (p + p_inf) / (1 / rho - b)) ** 0.5
    u = (p + gamma * p_inf) * (1 / rho - b) / (gamma - 1) + q
    h = u + p / rho
    state["u"] = u
    state["h"] = h
    state["Cv"] = Cv
    state["Cp"] = Cp
    state["gamma"] = gamma
    return state