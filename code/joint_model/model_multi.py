"""Multi-taxon joint density model (stage 2 of docs/JOINT_MODEL_DESIGN.md, section 4.5 test 2), in NumPyro.

EXPERIMENTAL, review only. Every taxon keeps the single-taxon model of model.py (same sites, same priors, same
likelihood; sites are prefixed "<tag>:"). The taxa are linked by lagged paths between latent states at the same station:

    latent log density of target k at station s, year t
        = [single-taxon terms for k]  +  sum over sources j of  beta_kj * c_j[s, t - lag]

    c_j[s, t] = (station level + station AR(1) deviation + protection function + paths into j) for source j,
                minus that station's own mean over its surveyed span.

So c_j is the source taxon's log density at that station relative to the station's usual level, with the region-year
effect of the source removed (regional swings such as a heatwave do not count as predictor variation). beta_kj is an
elasticity: the change in log density of k that follows a one-unit change in log density of j a year earlier.

Because the protection function of the source is part of c_j, a protection effect on predators carries through to
urchins, and from urchins to kelp. The response curve g_k of a target is then the DIRECT effect of protection (what is
not carried by the measured sources), and

    total_k(age) = g_k(age) + sum_j beta_kj * total_j(age - lag)

is the total. Both are reported (E10_direct, E10_total, E10_via:<source>).

Reading the paths as causal assumes the elasticity seen in station-level fluctuations and in differences between MPAs
also applies to the average change caused by protection, and that no unmeasured local driver moves both taxa.
"""
import numpy as np
import jax.numpy as jnp
from jax.scipy.special import gammaln
from jax import lax
import numpyro
import numpyro.distributions as dist
from numpyro.handlers import scope

HN = dist.HalfNormal
KNOTS = np.array([0, 4, 8, 12, 16, 20], dtype=float)


def basis_row(age):
    return np.clip((age - KNOTS[:-1]) / np.diff(KNOTS), 0.0, 1.0)


def latent(D, curve="rw", shock=True, centered=True, level_dev=False, marg_b=True, drift="regional", inc_scale=2.0,
           tau_d_scale=1.0, cluster_drift=False, station_trend=False, pretrend=False):
    """Single-taxon latent terms on the year x station grid. Mirrors model.py line for line (checked numerically by
    test_multi.py); returns level (station baseline + AR(1) deviation), region-year effects, protection function."""
    S, M, P, R, T = len(D["stations"]), len(D["mpas"]), len(D["programs"]), len(D["regions"]), len(D["years"])
    J = D["grid_B"].shape[-1]
    mu_b = numpyro.sample("mu_b", dist.Normal(D["log_mean_density"], 3.0).expand([P]))
    sigma_b = numpyro.sample("sigma_b", HN(1.5).expand([P]))
    marg_b = marg_b and centered
    if marg_b:
        b0 = mu_b[D["st_prog"]]
    else:
        b0 = mu_b[D["st_prog"]] + sigma_b[D["st_prog"]] * numpyro.sample("z_b", dist.Normal(0, 1).expand([S]))
    mu_lam = numpyro.sample("mu_lambda", dist.Normal(0, 1.0))
    sigma_lam = numpyro.sample("sigma_lambda", HN(1.0))
    z_lam = numpyro.sample("z_lambda", dist.Normal(0, 1).expand([M]))
    lam = jnp.where(jnp.asarray(D["mpa_ba"]) == 1, mu_lam + sigma_lam * z_lam, 3.0 * z_lam)
    numpyro.deterministic("lambda", lam)
    b = b0 + D["st_prot"] * lam[D["st_mpa"]]
    sigma_g = numpyro.sample("sigma_region_year", HN(1.0))
    gam = sigma_g * numpyro.sample("z_region_year", dist.Normal(0, 1).expand([R, T]))
    if drift != "none":
        s_a = numpyro.sample("sigma_drift", HN(0.3))
        alpha = jnp.cumsum(s_a * numpyro.sample("z_drift", dist.Normal(0, 1).expand([T])))
        alpha = alpha - alpha.mean()
        numpyro.deterministic("drift", alpha)
        gam = gam + alpha[None, :]
        if drift == "regional":
            s_r = numpyro.sample("sigma_drift_region", HN(0.3))
            beta = jnp.cumsum(s_r * numpyro.sample("z_drift_region", dist.Normal(0, 1).expand([R, T])), axis=1)
            gam = gam + beta - beta.mean(1, keepdims=True)
    rho = 0.95 * numpyro.sample("rho_raw", dist.Beta(2.0, 2.0))
    numpyro.deterministic("rho", rho)
    sigma_d = numpyro.sample("sigma_process", HN(1.0))
    first, inside = D["span_first"], D["span_inside"]
    if centered:
        delta = numpyro.sample("delta", dist.ImproperUniform(dist.constraints.real, (), (T, S)))
        lp_out = dist.Normal(0.0, 1.0).log_prob(delta)
        if marg_b:
            r = delta + D["log_mean_density"] - b[None, :]
            prev_r = jnp.concatenate([jnp.zeros((1, S)), r[:-1]], 0)
            lp_first = dist.Normal(0.0, sigma_d / jnp.sqrt(1 - rho ** 2)).log_prob(r)
            lp_next = dist.Normal(rho * prev_r, sigma_d).log_prob(r)
            has_next = jnp.concatenate([inside[1:], jnp.zeros((1, S))], 0)
            w = (first + inside) * ((1 - rho ** 2 * first) + rho ** 2 * has_next - rho * (inside + has_next)) / sigma_d ** 2
            a = jnp.sum(w, 0); c = jnp.sum(w * r, 0); sb2 = sigma_b[D["st_prog"]] ** 2
            rank_one = jnp.sum(-0.5 * jnp.log1p(sb2 * a) + 0.5 * sb2 * c ** 2 / (1 + sb2 * a))
            numpyro.factor("ar1", jnp.sum(first * lp_first + inside * lp_next + (1 - first - inside) * lp_out) + rank_one)
            level = delta + D["log_mean_density"]
        else:
            prev = jnp.concatenate([jnp.zeros((1, S)), delta[:-1]], 0)
            lp_first = dist.Normal(0.0, sigma_d / jnp.sqrt(1 - rho ** 2)).log_prob(delta)
            lp_next = dist.Normal(rho * prev, sigma_d).log_prob(delta)
            numpyro.factor("ar1", jnp.sum(first * lp_first + inside * lp_next + (1 - first - inside) * lp_out))
            level = b[None, :] + delta
    else:
        z = numpyro.sample("delta", dist.Normal(0, 1).expand([T, S]))
        def step(prev, x):
            zt, f, i = x
            cur = f * sigma_d * zt / jnp.sqrt(1 - rho ** 2) + i * (rho * prev + sigma_d * zt)
            return cur, cur
        _, delta = lax.scan(step, jnp.zeros(S), (z, first, inside))
        level = b[None, :] + delta
    if curve == "piecewise":
        s_inc = numpyro.sample("sigma_curve", HN(0.4 * inc_scale))
        inc = s_inc * numpyro.sample("z_curve", dist.Normal(0, 1).expand([J]))
    else:
        first_inc = numpyro.sample("inc_first", dist.Normal(0, 0.6 * inc_scale))
        s_inc = numpyro.sample("sigma_curve", HN(0.3))
        steps = s_inc * numpyro.sample("z_curve", dist.Normal(0, 1).expand([J - 1]))
        inc = first_inc + jnp.concatenate([jnp.zeros(1), jnp.cumsum(steps)])
    numpyro.deterministic("increments", inc)
    numpyro.deterministic("g_knots", jnp.cumsum(inc))
    tau_D = numpyro.sample("tau_D", HN(1.0 * tau_d_scale))
    Dm = tau_D * numpyro.sample("z_D", dist.Normal(0, 1).expand([M]))
    numpyro.deterministic("D_mpa", Dm)
    if level_dev:
        tau_L = numpyro.sample("tau_L", HN(1.0))
        Lm = tau_L * numpyro.sample("z_L", dist.Normal(0, 1).expand([M]))
        numpyro.deterministic("L_mpa", Lm)
    else:
        Lm = jnp.zeros(M)
    if shock:
        sigma_w = numpyro.sample("sigma_shock", HN(0.5))
        omega = sigma_w * numpyro.sample("z_shock", dist.Normal(0, 1).expand([T]))
        numpyro.deterministic("shock", omega)
    else:
        omega = jnp.zeros(T)
    m = D["st_mpa"]
    F = D["grid_after"] * (D["grid_B"] @ inc + Lm[m][None, :] + Dm[m][None, :] * D["grid_age"] / 10.0 + omega[:, None])
    # Structure options added 2026-11-11 after the placebo failures (see STRUCTURE_PLAN_v3.md). Both are off by default.
    local = jnp.zeros_like(F); trend = jnp.zeros_like(F)
    if cluster_drift:
        # Local trend: a random-walk year effect shared by one MPA's protected stations and the reference stations
        # assigned to it (design 4.3, the cluster-year term). The regional trend then only has to describe the region;
        # protection is judged against that MPA's own references.
        s_c = numpyro.sample("sigma_cluster", HN(0.3))
        chi = jnp.cumsum(s_c * numpyro.sample("z_cluster", dist.Normal(0, 1).expand([M, T])), axis=1)
        chi = chi - chi.mean(1, keepdims=True)
        local = chi[m].T
    if station_trend:
        # Station slope: every station may drift on its own (log units per decade), so that one station's long-run
        # drift is not read as a protection effect unless protected stations depart from the spread of all stations.
        s_k = numpyro.sample("sigma_trend", HN(0.5))
        kap = s_k * numpyro.sample("z_trend", dist.Normal(0, 1).expand([S]))
        yrs = jnp.asarray(D["years"], dtype=float)
        trend = kap[None, :] * ((yrs - yrs.mean()) / 10.0)[:, None]
    if pretrend:
        # Pre-trend sensitivity (the Nov 3 policy's alternative to the flat baseline): at MPAs with before data, the
        # protected stations may already be diverging from their references before designation (log units per decade).
        # That line is extended through the after years, so the protection function measures departure from it.
        mu_p = numpyro.sample("mu_pretrend", dist.Normal(0, 1.0))
        s_p = numpyro.sample("sigma_pretrend", HN(0.5))
        phi = mu_p + s_p * numpyro.sample("z_pretrend", dist.Normal(0, 1).expand([M]))
        numpyro.deterministic("pretrend", phi * jnp.asarray(D["mpa_ba"]))
        trend = trend + phi[m][None, :] * D["grid_rel"]
    return dict(level=level, gam=gam, F=F, inc=inc, local=local, trend=trend)


def observe(D, state):
    """Count likelihood for one taxon given its latent log density on the grid (region-year effect included)."""
    P, V = len(D["programs"]), len(D["visit_prog"])
    s, t = D["s"], D["t"]
    sigma_v = numpyro.sample("sigma_visit", HN(1.0).expand([P]))
    nu = sigma_v[D["visit_prog"]] * numpyro.sample("z_visit", dist.Normal(0, 1).expand([V]))
    zeta = jnp.concatenate([jnp.zeros(1), numpyro.sample("zone", dist.Normal(0, 1.0).expand([3]))])
    mon = numpyro.sample("month_raw", dist.Normal(0, 0.5).expand([12]))
    mon = mon - mon.mean()
    inv_sqrt_phi = numpyro.sample("inv_sqrt_phi", HN(1.0).expand([P]))
    phi = 1.0 / jnp.clip(inv_sqrt_phi, 1e-4, None) ** 2      # numerical guard, see model.py
    eta = state[t, s] + nu[D["v"]] + zeta[D["zone"]] + mon[D["month"]]
    mean = jnp.exp(eta + D["log_area"])
    ph = phi[D["obs_prog"]]
    lpm = jnp.log(ph + mean)
    log_nb0 = -ph * jnp.log1p(mean / ph)
    ll_cells = D["n_units"] * log_nb0 + D["cell_y"] * (eta + D["log_area"] - lpm)
    ll_disp = jnp.sum(D["uniq_mult"] * gammaln(D["uniq_y"] + phi[D["uniq_prog"]])) - jnp.sum(D["prog_units"] * gammaln(phi))
    numpyro.factor("y", jnp.sum(ll_cells) + ll_disp)


def mass_layer(M, curve="rw", mass_scale=1.0, fluct=None, mass_link=False, expose=False, mass_obs="gamma"):
    """Mean individual mass on the year x station grid of the stations that have measured animals (design 4.2).

        log M[s, t] = station level + siting difference (protected stations) + region-year effect (independent part,
                      shared random walk, regional random walk) + station AR(1) deviation
                      + after * (response curve by years since designation + MPA rate deviation * age / 10)
                      [+ theta * same-year fluctuation of the taxon's own density at that station, if mass_link]
        mean mass of a group of n animals ~ Gamma(n * k[program], mean M * zone effect)                 (mass_obs="gamma")
        or, mass_obs="lognormal": individual masses lognormal with log-SD s[program] and mean M * zone effect; the mean of
        n of them is taken as lognormal with the same mean and variance (exact for n = 1; Fenton-Wilkinson otherwise):
            s_n^2 = log(1 + (exp(s^2) - 1) / n),   log(group mean) ~ Normal(log(M * zone) - s_n^2 / 2, s_n)
        Added 2026-11-12 after the predictive check (M7): Gamma replicates of small groups have a much longer lower
        tail than the data. Off by default.

    Same response-curve machinery as density, with priors half as wide (a decade's change in mean mass beyond about
    x2.7 is implausible for slow-growing animals). Returns the response increments."""
    S, Mn, P, R, T = len(M["mass_stations"]), len(M["mpa_ba"]), int(M["st_prog"].max()) + 1, int(M["st_region"].max()) + 1, M["span"].shape[0]
    P = max(P, int(M["g_prog"].max()) + 1); J = M["grid_B"].shape[-1]
    mu_c = numpyro.sample("m_mu", dist.Normal(M["log_mean_mass"], 2.0).expand([P]))
    sigma_c = numpyro.sample("m_sigma_station", HN(1.0).expand([P]))
    c = mu_c[M["st_prog"]] + sigma_c[M["st_prog"]] * numpyro.sample("m_z_station", dist.Normal(0, 1).expand([S]))
    mu_lam = numpyro.sample("m_mu_lambda", dist.Normal(0, 0.5))
    sigma_lam = numpyro.sample("m_sigma_lambda", HN(0.5))
    z_lam = numpyro.sample("m_z_lambda", dist.Normal(0, 1).expand([Mn]))
    lam = jnp.where(jnp.asarray(M["mpa_ba"]) == 1, mu_lam + sigma_lam * z_lam, 1.5 * z_lam)
    c = c + M["st_prot"] * lam[M["st_mpa"]]
    gam = numpyro.sample("m_sigma_region_year", HN(0.5)) * numpyro.sample("m_z_region_year", dist.Normal(0, 1).expand([R, T]))
    alpha = jnp.cumsum(numpyro.sample("m_sigma_drift", HN(0.15)) * numpyro.sample("m_z_drift", dist.Normal(0, 1).expand([T])))
    beta = jnp.cumsum(numpyro.sample("m_sigma_drift_region", HN(0.15)) * numpyro.sample("m_z_drift_region", dist.Normal(0, 1).expand([R, T])), axis=1)
    gam = gam + (alpha - alpha.mean())[None, :] + beta - beta.mean(1, keepdims=True)
    rho = 0.95 * numpyro.sample("m_rho_raw", dist.Beta(2.0, 2.0))
    numpyro.deterministic("m_rho", rho)
    sigma_e = numpyro.sample("m_sigma_process", HN(0.5))
    z = numpyro.sample("m_delta", dist.Normal(0, 1).expand([T, S]))
    def step(prev, x):
        zt, f, i = x
        cur = f * sigma_e * zt / jnp.sqrt(1 - rho ** 2) + i * (rho * prev + sigma_e * zt)
        return cur, cur
    _, e = lax.scan(step, jnp.zeros(S), (z, M["span_first"], M["span_inside"]))
    if curve == "piecewise":
        s_inc = numpyro.sample("m_sigma_curve", HN(0.2 * mass_scale))
        inc = s_inc * numpyro.sample("m_z_curve", dist.Normal(0, 1).expand([J]))
    else:
        first_inc = numpyro.sample("m_inc_first", dist.Normal(0, 0.3 * mass_scale))
        s_inc = numpyro.sample("m_sigma_curve", HN(0.15 * mass_scale))
        inc = first_inc + jnp.concatenate([jnp.zeros(1), jnp.cumsum(s_inc * numpyro.sample("m_z_curve", dist.Normal(0, 1).expand([J - 1])))])
    numpyro.deterministic("m_increments", inc)
    numpyro.deterministic("m_g_knots", jnp.cumsum(inc))
    tau_D = numpyro.sample("m_tau_D", HN(0.5))
    Dm = tau_D * numpyro.sample("m_z_D", dist.Normal(0, 1).expand([Mn]))
    numpyro.deterministic("m_D_mpa", Dm)
    m = M["st_mpa"]
    logM = c[None, :] + gam[M["st_region"]].T + e + M["grid_after"] * (M["grid_B"] @ inc + Dm[m][None, :] * M["grid_age"] / 10.0)
    if mass_link and fluct is not None:
        theta = numpyro.sample("m_theta", dist.Normal(0, 0.5))
        logM = logM + theta * fluct[:, M["dens_index"]]
    zeta = jnp.concatenate([jnp.zeros(1), numpyro.sample("m_zone", dist.Normal(0, 0.5).expand([3]))])
    mean = jnp.exp(logM[M["g_t"], M["g_s"]] + zeta[M["g_zone"]])
    if mass_obs == "lognormal":
        sd = numpyro.sample("m_sd_log", HN(1.0).expand([P]))
        s2n = jnp.log1p(jnp.expm1(sd[M["g_prog"]] ** 2) / M["g_n"])
        numpyro.sample("m_obs", dist.LogNormal(jnp.log(mean) - 0.5 * s2n, jnp.sqrt(s2n)), obs=M["g_mean"])
    else:
        k = numpyro.sample("m_shape", dist.LogNormal(0.0, 1.0).expand([P]))
        conc = M["g_n"] * k[M["g_prog"]]
        numpyro.sample("m_obs", dist.Gamma(conc, conc / mean), obs=M["g_mean"])
    if expose:   # for predictive checks: expected mean mass of every fitted group and of every held-out group
        numpyro.deterministic("m_mu_groups", mean)
        if "h_s" in M and len(M["h_s"]):
            numpyro.deterministic("m_mu_holdout", jnp.exp(logM[M["h_t"], M["h_s"]] + zeta[M["h_zone"]]))
    return inc


def model(data, edges, config, curve="rw", lag=1, beta_scale=0.5, paths=True, split=False, mass=None, mass_scale=1.0, mass_link=False, expose=False, mass_obs="gamma", **common):
    """data: {tag: inputs}; edges: [(target, source, link)]; config: {tag: dict(centered=, level_dev=)}.
    paths=False fits the taxa side by side with no links (must reproduce the single-taxon fits).
    split=True gives every link two coefficients instead of one:
        beta      on the source's station-level fluctuation (baseline + AR(1) deviation + what its own sources'
                  fluctuations carry), centred within station. This is year-to-year co-variation at a station.
        beta_prot on the source's protection function (plus what its own sources' protection functions carry).
                  This is learned from how much the source responded to protection, MPA by MPA and year by year.
      Only beta_prot carries a protection effect downstream. With one coefficient (split=False) the two kinds of
      evidence are forced to agree; in the first real-data fit (2026-11-11) they did not, and the model bent the
      lobster protection estimate to reconcile them."""
    dev, prot, inc, total = {}, {}, {}, {}
    E = len(edges)
    if paths and E:
        beta = numpyro.sample("beta", dist.Normal(0, beta_scale).expand([E]))
        beta_p = numpyro.sample("beta_prot", dist.Normal(0, beta_scale).expand([E])) if split else beta
    for g, D in data.items():
        with scope(prefix=g, divider=":"):
            L = latent(D, curve=curve, **config[g], **common)
        path = jnp.zeros_like(L["F"]); path_p = jnp.zeros_like(L["F"]); via = {}
        if paths:
            for e, (k, j, lk) in enumerate(edges):
                if k == g:
                    path = path + beta[e] * lk["valid"] * dev[j][lk["src_t"]][:, lk["src_s"]]
                    if split:
                        path_p = path_p + beta_p[e] * lk["valid"] * prot[j][lk["src_t"]][:, lk["src_s"]]
                    via[j] = beta_p[e]
        n = jnp.sum(D["span"], 0)
        if split:
            prot[g] = L["F"] + path_p
            d = L["level"] + L["trend"] + path
            dev[g] = D["span"] * (d - jnp.sum(D["span"] * d, 0) / jnp.where(n > 0, n, 1.0))
            d = d + prot[g]
        else:
            d = L["level"] + L["trend"] + L["F"] + path
            dev[g] = D["span"] * (d - jnp.sum(D["span"] * d, 0) / jnp.where(n > 0, n, 1.0))
        inc[g] = L["inc"]
        # regional mean effect of protection at a given age: direct curve plus what the sources carry
        def tot(age, g=g, via=via):
            return jnp.asarray(basis_row(age)) @ inc[g] + sum(bj * total[j](max(age - lag, 0.0)) for j, bj in via.items())
        total[g] = tot
        with scope(prefix=g, divider=":"):
            numpyro.deterministic("E10_direct", jnp.asarray(basis_row(10.0)) @ L["inc"])
            numpyro.deterministic("E10_total", tot(10.0))
            for j, bj in via.items():
                numpyro.deterministic(f"E10_via_{j}", bj * total[j](10.0 - lag))
            observe(D, d + L["local"] + L["gam"][D["st_region"]].T)
            # mean individual mass and biomass (design 4.2): log biomass effect = log density effect + log mass effect
            if mass is not None and mass.get(g) is not None:
                lv = L["level"]
                fl = D["span"] * (lv - jnp.sum(D["span"] * lv, 0) / jnp.where(n > 0, n, 1.0)) if mass_link else None
                minc = mass_layer(mass[g], curve=curve, mass_scale=mass_scale, fluct=fl, mass_link=mass_link, expose=expose, mass_obs=mass_obs)
                e10m = jnp.asarray(basis_row(10.0)) @ minc
                numpyro.deterministic("E10_mass", e10m)
                numpyro.deterministic("E10_biomass", tot(10.0) + e10m)
