"""Checks of the multi-taxon model code.
1. With no links, each taxon's part of the multi-taxon log density equals the single-taxon model's (model.py)
   at the same parameter values, for both sampler parameterisations.
2. The link maps point at the same station, `lag` years earlier, inside the source's surveyed span.
Run: python3 code/joint_model/test_multi.py [path to supported_measurement_uses.parquet]
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import run as single
from run import attach_counts, USES
import jax
from numpyro.infer.util import initialize_model
from numpyro.infer import init_to_uniform
from data import age_basis
from data_multi import build_multi, CONFIG, TAXA
import model as m1
import model_multi as mm

uses = sys.argv[1] if len(sys.argv) > 1 else USES
tags = ["lob", "shp", "pur", "red", "kelp"]
data, edges = build_multi(uses, tags)
worst = 0.0
for g in tags:
    D = attach_counts(dict(data[g]), data[g]["y"]); D.pop("cells")
    B = age_basis(D["age"])
    for curve in ("rw", "piecewise"):
        kw1 = dict(curve=curve, shock=True, centered=CONFIG[g]["centered"], level_dev=CONFIG[g]["level_dev"], marg_b=True,
                   drift="regional", inc_scale=2.0)
        i1 = initialize_model(jax.random.PRNGKey(3), m1.model, init_strategy=init_to_uniform(radius=0.5),
                              model_args=(D, B), model_kwargs=dict(y=D["y"], **kw1))
        i2 = initialize_model(jax.random.PRNGKey(3), mm.model, init_strategy=init_to_uniform(radius=0.5),
                              model_args=({g: D}, [], {g: CONFIG[g]}), model_kwargs=dict(curve=curve, paths=False, drift="regional", inc_scale=2.0))
        rng = np.random.default_rng(7)
        z1 = {k: 0.5 * rng.normal(size=np.shape(v)) for k, v in i1.param_info.z.items()}
        z2 = {f"{g}:{k}": v for k, v in z1.items()}
        assert set(z2) == set(i2.param_info.z), (set(z2) ^ set(i2.param_info.z))
        p1, p2 = float(i1.potential_fn(z1)), float(i2.potential_fn(z2))
        worst = max(worst, abs(p1 - p2) / abs(p1))
        print(f"{g:5s} {curve:9s} potential single {p1:.6f} multi {p2:.6f} diff {p1 - p2:.2e}")
        assert abs(p1 - p2) < 1e-6 * abs(p1)
print("1. log density matches the single-taxon model; worst relative difference", worst)

for k, j, lk in edges:
    Dk, Dj = data[k], data[j]
    rng = np.random.default_rng(1); n = 0
    for _ in range(4000):
        t, s = rng.integers(len(Dk["years"])), rng.integers(len(Dk["stations"]))
        if lk["valid"][t, s]:
            assert Dj["stations"][lk["src_s"][s]] == Dk["stations"][s]
            assert Dj["years"][lk["src_t"][t]] == Dk["years"][t] - 1
            assert Dj["span"][lk["src_t"][t], lk["src_s"][s]] == 1
            n += 1
    print(f"2. {k}<-{j}: {lk['n_linked']} of {lk['n_target_station_years']} target station-years linked, "
          f"{lk['n_stations_linked']} stations; {n} sampled links checked")
print("ALL CHECKS PASSED")

# 3. Mass layer: its contribution to the joint log density (priors + likelihood) equals an independent
#    NumPy / SciPy calculation at the same parameter values, for both curve priors.
from scipy import stats as st
from numpyro.infer.util import log_density, constrain_fn
from data_mass import build_mass
from model_multi import basis_row


def mass_logp(M, P_, curve, mass_scale=1.0, obs="gamma"):
    v = lambda k: np.asarray(P_[k], dtype=float)
    S, T = len(M["mass_stations"]), M["span"].shape[0]
    lp = st.norm(M["log_mean_mass"], 2.0).logpdf(v("m_mu")).sum() + st.halfnorm(scale=1.0).logpdf(v("m_sigma_station")).sum()
    lp += st.norm.logpdf(v("m_z_station")).sum() + st.norm(0, 0.5).logpdf(v("m_mu_lambda")) + st.halfnorm(scale=0.5).logpdf(v("m_sigma_lambda"))
    lp += st.norm.logpdf(v("m_z_lambda")).sum() + st.halfnorm(scale=0.5).logpdf(v("m_sigma_region_year")) + st.norm.logpdf(v("m_z_region_year")).sum()
    lp += st.halfnorm(scale=0.15).logpdf(v("m_sigma_drift")) + st.norm.logpdf(v("m_z_drift")).sum()
    lp += st.halfnorm(scale=0.15).logpdf(v("m_sigma_drift_region")) + st.norm.logpdf(v("m_z_drift_region")).sum()
    lp += st.beta(2, 2).logpdf(v("m_rho_raw")) + st.halfnorm(scale=0.5).logpdf(v("m_sigma_process")) + st.norm.logpdf(v("m_delta")).sum()
    if curve == "piecewise":
        lp += st.halfnorm(scale=0.2 * mass_scale).logpdf(v("m_sigma_curve")) + st.norm.logpdf(v("m_z_curve")).sum()
        inc = v("m_sigma_curve") * v("m_z_curve")
    else:
        lp += st.norm(0, 0.3 * mass_scale).logpdf(v("m_inc_first")) + st.halfnorm(scale=0.15 * mass_scale).logpdf(v("m_sigma_curve")) + st.norm.logpdf(v("m_z_curve")).sum()
        inc = v("m_inc_first") + np.concatenate([[0.0], np.cumsum(v("m_sigma_curve") * v("m_z_curve"))])
    lp += st.halfnorm(scale=0.5).logpdf(v("m_tau_D")) + st.norm.logpdf(v("m_z_D")).sum() + st.norm(0, 0.5).logpdf(v("m_zone")).sum()
    lp += (st.halfnorm(scale=1.0).logpdf(v("m_sd_log")).sum() if obs == "lognormal" else st.lognorm(s=1.0, scale=1.0).logpdf(v("m_shape")).sum())
    c = v("m_mu")[M["st_prog"]] + v("m_sigma_station")[M["st_prog"]] * v("m_z_station")
    lam = np.where(M["mpa_ba"] == 1, v("m_mu_lambda") + v("m_sigma_lambda") * v("m_z_lambda"), 1.5 * v("m_z_lambda"))
    c = c + M["st_prot"] * lam[M["st_mpa"]]
    al = np.cumsum(v("m_sigma_drift") * v("m_z_drift")); br = np.cumsum(v("m_sigma_drift_region") * v("m_z_drift_region"), axis=1)
    gam = v("m_sigma_region_year") * v("m_z_region_year") + (al - al.mean())[None, :] + br - br.mean(1, keepdims=True)
    rho, sig, z = 0.95 * v("m_rho_raw"), v("m_sigma_process"), v("m_delta")
    e = np.zeros((T, S))
    for s_ in range(S):                       # AR(1) over each station's own span, written out as a loop
        prev = 0.0
        for t_ in range(T):
            if M["span_first"][t_, s_]:
                prev = sig * z[t_, s_] / np.sqrt(1 - rho ** 2)
            elif M["span_inside"][t_, s_]:
                prev = rho * prev + sig * z[t_, s_]
            else:
                prev = 0.0
            e[t_, s_] = prev
    Dm = v("m_tau_D") * v("m_z_D")
    zeta = np.concatenate([[0.0], v("m_zone")])
    lik = 0.0
    for i in range(len(M["g_n"])):            # one group at a time
        t_, s_ = M["g_t"][i], M["g_s"][i]
        age = M["grid_age"][t_, s_]
        logM = c[s_] + gam[M["st_region"][s_], t_] + e[t_, s_] + M["grid_after"][t_, s_] * (basis_row(age) @ inc + Dm[M["st_mpa"][s_]] * age / 10.0)
        mean = np.exp(logM + zeta[M["g_zone"][i]])
        if obs == "lognormal":                # mean of n lognormal animals, matched on mean and variance
            var_one = np.exp(v("m_sd_log")[M["g_prog"][i]] ** 2) - 1.0
            s2 = np.log(1.0 + var_one / M["g_n"][i])
            lik += st.lognorm(s=np.sqrt(s2), scale=mean * np.exp(-0.5 * s2)).logpdf(M["g_mean"][i])
        else:
            k = M["g_n"][i] * v("m_shape")[M["g_prog"][i]]
            lik += st.gamma(a=k, scale=mean / k).logpdf(M["g_mean"][i])
    return float(lp + lik)


worst = 0.0
for g in ("lob", "pur", "shp"):
    D = attach_counts(dict(data[g]), data[g]["y"]); D.pop("cells")
    M = build_mass(uses, data[g], TAXA[g])
    for curve, obs in (("rw", "gamma"), ("piecewise", "gamma"), ("rw", "lognormal")):
        args = ({g: D}, [], {g: CONFIG[g]}); kw = dict(curve=curve, paths=False, drift="regional", inc_scale=2.0)
        mk = dict(mass={g: M}, mass_obs=obs, **kw)
        i2 = initialize_model(jax.random.PRNGKey(5), mm.model, init_strategy=init_to_uniform(radius=0.5), model_args=args, model_kwargs=mk)
        rng = np.random.default_rng(11)
        z = {k: 0.3 * rng.normal(size=np.shape(v_)) for k, v_ in sorted(i2.param_info.z.items())}
        con = constrain_fn(mm.model, args, mk, z)
        with_mass = float(log_density(mm.model, args, mk, con)[0])
        without = float(log_density(mm.model, args, kw, {k: v_ for k, v_ in con.items() if ":m_" not in k})[0])
        ind = mass_logp(M, {k.split(":")[1]: v_ for k, v_ in con.items() if ":m_" in k}, curve, obs=obs)
        d = abs((with_mass - without) - ind) / abs(ind); worst = max(worst, d)
        print(f"3. {g:4s} {curve:9s} {obs:9s} mass layer log density: model {with_mass - without:.6f}  independent {ind:.6f}  ({M['n_groups']} groups)")
        assert d < 1e-9, d
print("3. mass layer matches the independent calculation; worst relative difference", worst)

# 4. Canopy (Landsat) likelihood: the model's observation term at random parameter values equals a separate SciPy
#    calculation, one polygon-quarter at a time. Skipped with a note when the canopy table is not available.
from data_canopy import build_canopy, find_canopy, attach_canopy
from numpyro.handlers import substitute, trace
cpath = find_canopy(uses)
if cpath is None:
    print("4. canopy table not found: canopy check SKIPPED")
    print("ALL CHECKS PASSED (1-3)")
else:
    from data_multi import grid
    for season in ("spring", "all", "max"):
        Dc = grid(build_canopy(cpath, season=season)); Dc.pop("cells"); Dc = attach_canopy(Dc, Dc["y"])
        args = ({"can": Dc}, [], {"can": CONFIG["can"]}); kw = dict(curve="rw", paths=False, drift="regional", inc_scale=2.0)
        i4 = initialize_model(jax.random.PRNGKey(7), mm.model, init_strategy=init_to_uniform(radius=0.5), model_args=args, model_kwargs=kw)
        rng = np.random.default_rng(13)
        con = constrain_fn(mm.model, args, kw, {k: 0.3 * rng.normal(size=np.shape(v_)) for k, v_ in sorted(i4.param_info.z.items())})
        full = float(log_density(mm.model, args, kw, con)[0])
        empty = dict(Dc, c_s=Dc["c_s"][:0], c_t=Dc["c_t"][:0], c_q=Dc["c_q"][:0], c_logpix=Dc["c_logpix"][:0], c_pos=Dc["c_pos"][:0], c_logy=Dc["c_logy"][:0])
        none = float(log_density(mm.model, ({"can": empty}, [], {"can": CONFIG["can"]}), kw, con)[0])
        mu = np.asarray(trace(substitute(mm.model, data=con)).get_trace(*args, **dict(kw, expose=True))["can:c_mu"]["value"], dtype=float)
        a_, b_, sd_ = float(con["can:c_det0"]), float(con["can:c_det_slope"]), float(con["can:c_sd"])
        ind = 0.0
        for i in range(len(mu)):
            p_det = 1.0 / (1.0 + np.exp(-(a_ + b_ * (mu[i] + Dc["c_logpix"][i] - Dc["c_centre"]))))
            ind += (np.log(p_det) + st.t(df=5, loc=mu[i], scale=sd_).logpdf(Dc["c_logy"][i])) if Dc["c_pos"][i] else np.log1p(-p_det)
        d = abs((full - none) - ind) / abs(ind)
        print(f"4. canopy {season:6s} observation term: model {full - none:.6f}  independent {ind:.6f}  ({len(mu)} polygon-quarters, {int((1 - Dc['c_pos']).sum())} without canopy)")
        assert d < 1e-9, d
    print("ALL CHECKS PASSED (1-4)")
