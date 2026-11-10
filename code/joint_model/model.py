"""Joint hierarchical density model (stage 1 of docs/JOINT_MODEL_DESIGN.md), in NumPyro.

Counts per visit x zone ~ negative binomial with area offset. Latent log density per station-year =
station baseline + region-year effect + station AR(1) + protection function (protected stations only).
Protection function = siting difference + after * (regional response curve by years since designation
+ MPA-specific rate deviation + protected-by-year shock shared across MPAs).
"""
import jax.numpy as jnp
from jax.scipy.special import gammaln
from jax import lax
import numpyro
import numpyro.distributions as dist


def model(D, B, curve="piecewise", shock=True, y=None, centered=True, level_dev=False):
    S, M, P, R, T = len(D["stations"]), len(D["mpas"]), len(D["programs"]), len(D["regions"]), len(D["years"])
    V = len(D["visit_prog"]); J = B.shape[1]
    HN = dist.HalfNormal
    # station baselines
    mu_b = numpyro.sample("mu_b", dist.Normal(D["log_mean_density"], 3.0).expand([P]))
    sigma_b = numpyro.sample("sigma_b", HN(1.5).expand([P]))
    z_b = numpyro.sample("z_b", dist.Normal(0, 1).expand([S]))
    # siting difference
    mu_lam = numpyro.sample("mu_lambda", dist.Normal(0, 1.0))
    sigma_lam = numpyro.sample("sigma_lambda", HN(1.0))
    z_lam = numpyro.sample("z_lambda", dist.Normal(0, 1).expand([M]))
    ba = jnp.asarray(D["mpa_ba"])
    lam = jnp.where(ba == 1, mu_lam + sigma_lam * z_lam, 3.0 * z_lam)
    numpyro.deterministic("lambda", lam)
    b = mu_b[D["st_prog"]] + sigma_b[D["st_prog"]] * z_b + D["st_prot"] * lam[D["st_mpa"]]
    # region-year
    sigma_g = numpyro.sample("sigma_region_year", HN(1.0))
    gam = sigma_g * numpyro.sample("z_region_year", dist.Normal(0, 1).expand([R, T]))
    # station AR(1): states sampled directly over each station's observed span (centered parameterization)
    rho = 0.95 * numpyro.sample("rho_raw", dist.Beta(2.0, 2.0))
    numpyro.deterministic("rho", rho)
    sigma_d = numpyro.sample("sigma_process", HN(1.0))
    first, inside = D["span_first"], D["span_inside"]            # T x S masks
    if centered:   # best when station-years are well informed and the process SD is clearly above zero (urchins)
        delta = numpyro.sample("delta", dist.ImproperUniform(dist.constraints.real, (), (T, S)))
        prev = jnp.concatenate([jnp.zeros((1, S)), delta[:-1]], 0)
        lp_first = dist.Normal(0.0, sigma_d / jnp.sqrt(1 - rho ** 2)).log_prob(delta)
        lp_next = dist.Normal(rho * prev, sigma_d).log_prob(delta)
        lp_out = dist.Normal(0.0, 1.0).log_prob(delta)
        numpyro.factor("ar1", jnp.sum(first * lp_first + inside * lp_next + (1 - first - inside) * lp_out))
    else:          # innovations parameterization: needed when the process SD is near zero (lobster, sheephead)
        z = numpyro.sample("delta", dist.Normal(0, 1).expand([T, S]))
        def step(prev, x):
            zt, f, i = x
            cur = f * sigma_d * zt / jnp.sqrt(1 - rho ** 2) + i * (rho * prev + sigma_d * zt)
            return cur, cur
        _, delta = lax.scan(step, jnp.zeros(S), (z, first, inside))
    # response
    if curve in ("piecewise", "rw"):
        if curve == "piecewise":   # independent increments: shrinks the curve toward zero
            s_inc = numpyro.sample("sigma_curve", HN(0.4))
            inc = s_inc * numpyro.sample("z_curve", dist.Normal(0, 1).expand([J]))
        else:                      # random-walk increments: shrinks the curve toward a straight line (design v2)
            first_inc = numpyro.sample("inc_first", dist.Normal(0, 0.6))
            s_inc = numpyro.sample("sigma_curve", HN(0.3))
            steps = s_inc * numpyro.sample("z_curve", dist.Normal(0, 1).expand([J - 1]))
            inc = first_inc + jnp.concatenate([jnp.zeros(1), jnp.cumsum(steps)])
        numpyro.deterministic("increments", inc)
        g_obs = B @ inc
        numpyro.deterministic("g_knots", jnp.cumsum(inc))
        numpyro.deterministic("E10", inc[0] + inc[1] + 0.5 * inc[2])
    else:
        e10 = numpyro.sample("E10", dist.Normal(0, 1.5))
        g_obs = e10 * D["age"] / 10.0
    tau_D = numpyro.sample("tau_D", HN(1.0))
    Dm = tau_D * numpyro.sample("z_D", dist.Normal(0, 1).expand([M]))
    numpyro.deterministic("D_mpa", Dm)
    if level_dev:  # MPA-specific shift in the level of the response (needed when the response is step-like)
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
    s, t = D["s"], D["t"]
    m_obs = D["st_mpa"][s]
    F = D["after"] * (g_obs + Lm[m_obs] + Dm[m_obs] * D["age"] / 10.0 + omega[t])
    # observation level
    sigma_v = numpyro.sample("sigma_visit", HN(1.0).expand([P]))
    nu = sigma_v[D["visit_prog"]] * numpyro.sample("z_visit", dist.Normal(0, 1).expand([V]))
    zeta = jnp.concatenate([jnp.zeros(1), numpyro.sample("zone", dist.Normal(0, 1.0).expand([3]))])
    mon = numpyro.sample("month_raw", dist.Normal(0, 0.5).expand([12]))
    mon = mon - mon.mean()
    inv_sqrt_phi = numpyro.sample("inv_sqrt_phi", HN(1.0).expand([P]))
    phi = 1.0 / inv_sqrt_phi ** 2
    eta = b[s] + gam[D["st_region"][s], t] + delta[t, s] + F + nu[D["v"]] + zeta[D["zone"]] + mon[D["month"]]
    mean = jnp.exp(eta + D["log_area"])
    # Negative binomial log likelihood, exact, in two parts: cell totals carry the mean; unit counts carry the
    # dispersion (lgamma evaluated once per distinct count and program). Constant lgamma(y + 1) is dropped.
    ph = phi[D["obs_prog"]]
    lpm = jnp.log(ph + mean)
    ll_cells = D["n_units"] * ph * (jnp.log(ph) - lpm) + D["cell_y"] * (eta + D["log_area"] - lpm)
    ll_disp = jnp.sum(D["uniq_mult"] * gammaln(D["uniq_y"] + phi[D["uniq_prog"]])) - jnp.sum(D["prog_units"] * gammaln(phi))
    numpyro.factor("y", jnp.sum(ll_cells) + ll_disp)
