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
