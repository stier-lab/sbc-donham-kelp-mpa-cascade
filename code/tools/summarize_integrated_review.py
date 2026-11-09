"""Plain-language review notes derived from the completed integrated tables."""
from pathlib import Path
import pandas as pd
import numpy as np

def write_summary():
    out=Path('outputs/integrated_review_20261103')
    p=pd.read_csv(out/'pooled_results.csv')
    p=p[p.Source.eq('All programs')]
    def result(taxon,response,model,term):
        r=p[(p.Taxa==taxon)&(p.Resp==response)&(p.model==model)&(p.term==term)].iloc[0]
        a,b,c=100*np.expm1([r.estimate,r.low,r.high])
        return f'{a:+.1f}% (95% interval {b:+.1f}% to {c:+.1f}%)'
    mc=pd.read_csv(out/'monte_carlo_check.csv')
    gained=pd.read_csv(out/'seven_day_recovered_years.csv')
    text=f'''# Reading the integrated review

This is a review of the corrected data and frozen rules. The manuscript and supplement have not changed.

## Coverage and the seven-day decision

The candidate inventory retains all 31 researcher/ecological exclusions. Of 175 retained histories, 62 support before–after models, another 100 support post-protection trends only, and 13 lack enough supported post-protection years. Before–after histories also have separately estimated post-protection trends. Missing sizes do not erase valid density observations.

Allowing same-site/year KFM urchin and kelp sizes within seven days recovers {len(gained)} comparison-years relative to exact-date matching. These are additional observations within histories, not {len(gained)} new sites. See [the recovered years](seven_day_recovered_years.csv) and [all fitted sensitivities](pooled_coverage_sensitivities.csv).

## What the patterns say

- After protection, the kelp protected/reference ratio increases by {result('M. pyrifera','Bio','post_trend','recovery_slope')} per year.
- Purple-urchin density becomes lower relative to reference sites: {result('S. purpuratus','Den','post_trend','recovery_slope')} per year.
- Red-urchin density also has a negative estimated trend, but its interval endpoint is close enough to zero to change sides between random seeds. Do not turn that small numerical distinction into a different biological conclusion.
- The pooled sheephead and lobster post-protection intervals include both increases and decreases. Individual programs and histories can differ; use the source-specific results rather than interpreting the pooled interval as no response anywhere.

These trends describe change relative to reference sites. They do not establish that protection caused those changes or that a predator-driven cascade occurred.

## Before–after interpretation

The flat-baseline model estimates a kelp ratio change of {result('M. pyrifera','Bio','flat_before_after','average_observed_post_change')} over observed post-protection years. Allowing the before-protection trend to continue gives {result('M. pyrifera','Bio','pretrend_sensitivity','average_observed_post_change')}. That difference matters: the inference depends on the assumed trajectory without protection.

The corresponding flat-baseline lobster density estimate is {result('P. interruptus','Den','flat_before_after','average_observed_post_change')}. Its positive estimate is imprecise. Before–after estimates and post-protection trends answer different questions and remain separate.

## Numerical and release checks

Across 134 pooled summaries, changing the measurement-draw seed moves point estimates by at most {mc.point_change_in_main_SE.max():.2f} standard errors. The earlier large KFM lobster instability was traced to finite-draw correlations between independent survey years. Known-zero measurement covariances are now enforced; biological temporal dependence is still estimated. Variance fitting uses multiple starts and a checked temporal-correlation profile.

Intervals remain conditional on the represented sampling units, conversion equations and fitted covariance. Group-mean sizes, allometry, detection, sparse/zero samples and missing habitat information limit what the files can establish. Boundary flags remain in the diagnostics. The [annual input decomposition](input_change_decomposition.csv) separates coverage and measurement changes using the same simple diagnostic trend; those slopes are not substitutes for the joint model results.

Before updating the paper, review these results and the remaining offset, satellite-timing and habitat limitations, then regenerate downstream pathway/resilience analyses, figures and concordance checks. KFM redistribution permission remains a separate requirement. The current tables do not certify a final paper release.
'''
    (out/'REVIEW_NOTES.md').write_text(text)
