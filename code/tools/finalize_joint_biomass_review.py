"""Verify review outputs, summarize sensitivity, then refresh report and guide."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from report_joint_borrowed_biomass import run as report
O=Path('outputs/joint_borrowed_biomass_20261102')
def run():
    base=pd.read_csv(O/'joint_trends.csv');b=base[base.rho.eq(0)].set_index('series_id');sens=[]
    for name in ['seed2','pool10','strict0']:
        assert (O/name/'completed.json').exists()
        for r in pd.read_csv(O/name/'joint_trends.csv').itertuples():
            ref=b.loc[r.series_id];sens.append(dict(series_id=r.series_id,scenario=name,n_years=r.n_years,estimate=r.estimate,low=r.low,high=r.high,estimate_change=r.estimate-ref.estimate,low_change=r.low-ref.low,high_change=r.high-ref.high))
    pd.DataFrame(sens).to_csv(O/'sensitivity_comparison.csv',index=False)
    for folder in [O,*[O/x for x in ['seed2','pool10','strict0']]]:
        hashes=json.loads((folder/'input_hashes.json').read_text())
        for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
        fits=pd.read_csv(folder/'bootstrap_fits.csv');assert fits.converged.all()
    assert len(base)==18 and len(pd.read_csv(O/'support.csv'))==6 and pd.read_csv(O/'support.csv').n_years.sum()==45
    ident=pd.read_csv(O/'identification.csv');assert (ident.active_components==ident.projected_kernel_rank).all()
    assert pd.read_csv(O/'multiple_starts.csv').slope_difference.abs().max()<1e-4
    tests=json.loads((O/'independent_checks.json').read_text());assert tests['success']
    assert json.loads((O/'diagnostic_checks.json').read_text())['execution_complete']
    report()
    summary=pd.read_csv(O/'simulation_summary.csv');warning=bool((summary.coverage_mc_high<.95).any())
    record=dict(review_outputs_verified=True,nominal_interval_coverage_certified=False,conditional_simulation_warning=warning,primary_updated=False,raw_support_unchanged=True,model_source_sha256=hashlib.sha256(Path('code/tools/joint_borrowed_biomass.py').read_bytes()).hexdigest())
    record['output_hashes']={name:hashlib.sha256((O/name).read_bytes()).hexdigest() for name in ['joint_trends.csv','same_year_comparison.csv','diagnostics_summary.csv','sensitivity_comparison.csv','simulation_summary.csv','independent_checks.json','diagnostic_checks.json','support.csv']}
    (O/'report_verified.json').write_text(json.dumps(record,indent=2))
    older=Path('outputs/biomass_next_steps_20261102/README.md')
    note='> Earlier review stage. The [joint borrowed-biomass model review](../joint_borrowed_biomass_20261102/README.md) now extends the uncertainty work below; retain this report as dated evidence.\n\n'
    if not older.read_text().startswith('> Earlier review stage.'):
        older.write_text(note+older.read_text())
    print(json.dumps(record))
if __name__=='__main__':run()
