"""Focused current-decision viewer; historical audits stay in their original files."""
from pathlib import Path
import json
import pandas as pd

def build(payload):
    summaries={
      'R01':'Keep sites excluded when researchers or ecological criteria call for it. Do not choose sites by the results.',
      'R02':'Use at least three matched, usable years before designation and three usable years after. Count the designation year as after. Keep all supported years, with no nine-year limit. Check four before years and a continuing pre-protection trend as alternatives.',
      'R03':'Compare protected and reference surveys only where they share recorded depth zones. Also check results using the same zones throughout.',
      'R04':'Keep valid surveys from any month. Give repeat visits equal weight within zones, then average zones within stations and stations within each side.',
      'R05':'Keep the 2004 CRANE fish surveys and average repeat surveys. Leave out incomplete passes when calculating the combined fish measure.',
      'R06':'Use the documented stations surveyed each year. Check whether results change when the station set stays constant.',
      'R07':'Treat missing counts and sizes as missing. Keep a density result when counts are sound but biomass is not.',
      'R08':'Keep excluded sites out, even if their names or group labels differ across files.',
      'R09':'Use sizes from the counted site and survey when possible. For KFM, allow sizes measured within seven days; use the separate tested borrowing rule for PISCO.',
      'R10':'Borrow sizes only from reviewed reference surveys of the same species, protected area, zone and year. Test each borrowing rule without using the target site to calibrate itself.',
      'R11':'For PISCO borrowing, allow no larger count–size date gap than the one tested for that comparison. Require at least two calibration site-years.',
      'R12':'Require biomass support for every retained station and visit in a comparison-year. Do not drop a station only to make biomass calculable.',
      'R13':'Account for repeated use of counts, sizes and reference sites, and for change through time. State what the available data cannot measure.',
      'R14':'Review these results before changing the paper. Then rerun the linked analyses and check every reported result.',
      'R15':'Set borrowing rules before looking at the effects. Show which biomass results depend on borrowed sizes and what happens with direct sizes only.'}
    reasons={
      'R01':'Researcher advice and ecological fit determine site eligibility, not whether a result is large or statistically clear.',
      'R02':'The same time rule should apply to every comparison. The alternative fits show how much the baseline choice matters.',
      'R03':'A depth difference can look like a biological change. A shared zone name helps, but it does not guarantee the same depth or habitat.',
      'R04':'There is no supported reason to discard all surveys outside one season. Repeated visits should not count as extra stations.',
      'R05':'A missing fish survey level is not a count of zero. The available CRANE records follow the approved survey rule.',
      'R06':'A changing set of stations could change the apparent trend.',
      'R07':'Missing size data should not erase a valid count history.',
      'R08':'Different names for the same excluded site must not bring its measurements back into the analysis.',
      'R09':'Sizes from other places or distant dates may not represent the animals counted.',
      'R10':'Borrowed sizes can alter biomass; donors should be comparable to the count survey.',
      'R11':'A larger date gap would go beyond what the data used to test borrowing can support. Two calibration site-years remain sparse evidence.',
      'R12':'Dropping unsupported stations only for biomass would change the site comparison.',
      'R13':'Shared measurements and reference sites make estimates less independent than they first appear.',
      'R14':'The paper should report one checked set of analyses and results.',
      'R15':'The direct-size-only check shows how much the biomass result relies on borrowing.'}
    rules=[dict(decision_id=r['rule_id'],decision=summaries[r['rule_id']],why=reasons[r['rule_id']],applies_to=r['scope'],full_rule=r['current_rule'],status=r['implementation_status'],limitation=r['limitation_or_supersedes']) for r in payload['Current rules']]
    pd.DataFrame(rules).to_csv('outputs/dataset_inventory_20261101/current_decisions_and_reasons.csv',index=False)
    select=lambda rows,keys:[{k:r.get(k) for k in keys} for r in rows]
    sections={
      'Decisions and reasons':rules,
      'Included and excluded':select(payload['Eligibility decisions'],['series_id','exclude','reason','mapped_protected_sites','mapped_reference_sites','reference_evidence','integrated_before_years','integrated_after_years','integrated_supported_years','integrated_eligible_model','reference_sensitivity']),
      'Survey coverage':payload['Integrated annual roster'],
      'Biomass calculations':payload['Biomass calculations'],
      'Biomass support':payload['Integrated biomass visits'],
      'Biomass coverage':json.loads(pd.read_csv('outputs/biomass_policy_review_20261103/comparison_support.csv').to_json(orient='records')),
      'Biomass sensitivities':json.loads(pd.read_csv('outputs/biomass_policy_review_20261103/pooled_sensitivity_results.csv').to_json(orient='records')),
      'Time series':payload['Integrated model eligibility'],
      'Updated results':json.loads(pd.read_csv('outputs/whole_dataset_biomass_update_20261105/pooled_effects.csv').to_json(orient='records')),
      'Change steps':json.loads(pd.read_csv('outputs/whole_dataset_biomass_update_20261105/series_change_steps.csv').to_json(orient='records')),
      'Pathway refits':json.loads(pd.read_csv('outputs/whole_dataset_biomass_update_20261105/matched_year_pathway_fits.csv').to_json(orient='records')),
      'Heatwave refits':json.loads(pd.read_csv('outputs/whole_dataset_biomass_update_20261105/heatwave_period_refits.csv').to_json(orient='records')),
      'Open questions':[r for r in payload['Questions'] if r['issue_id'] not in ['U06','U08','U14','U17']],
      'Uncertainty limits':payload['Integrated uncertainty limits']}
    meta={
      'Decisions and reasons':['Start here','The rules we are using now and why. Select “Show all columns” for the exact technical rule and its limits. These have not yet been put in the manuscript.',['decision_id','decision','why','applies_to']],
      'Included and excluded':['Sites and surveys','Each row is one program, protected–reference comparison, species and measurement. See why it is included or excluded, and which stations are involved. An included comparison may still lack enough years for a model.',['series_id','exclude','reason','mapped_protected_sites','mapped_reference_sites','integrated_eligible_model']],
      'Survey coverage':['Sites and surveys','Each row shows one comparison in one year: which stations were surveyed, whether both sides share a depth zone, and whether that year can be used.',['series_id','year','survey_eligible','reason','protected_stations','reference_stations','shared_zones']],
      'Biomass calculations':['Biomass','How we turn counts and sizes into biomass, including grouped sizes, missing lengths and actual zero counts.',['input','calculation','formula','scope','limitation']],
      'Biomass support':['Biomass','Each row shows a survey record used in a biomass check. One size sample can be used in several comparisons. When sizes are missing, the count may still support density.',['series_id','site','year','date','biomass_supported','basis','gap_days']],
      'Biomass coverage':['Biomass','For each comparison, see how many surveyed years have enough information for biomass, how many use borrowed sizes, and how many remain unsupported.',['series_id','candidate_years','supported_years','supported_borrowed_years','unsupported_years']],
      'Biomass sensitivities':['Biomass','Compare trends using the current borrowing rule, direct sizes only, and an assumed range for LTER group sizes. Direct sizes alone cover fewer years; the intervals do not include every source of uncertainty.',['Taxa','scenario','estimate_pct_per_year','low_pct_per_year','high_pct_per_year','comparisons']],
      'Time series':['Time series','See the years available before and after designation and which analysis each comparison supports. Three matched earlier years allow a before–after comparison; with fewer, we report only the later trend if at least three later years are available. We use all observed years. These patterns alone do not prove that protection caused the change.',['series_id','before_years','after_years','first_year','last_year','eligible_model','missing_response_years']],
      'Updated results':['Results for review','Pooled estimates from the corrected whole-dataset run. Post-protection slopes and before–after changes are separate; these results have not replaced the paper.',['Taxa','Resp','estimand','comparisons','estimate_percent','low_percent','high_percent']],
      'Change steps':['Results for review','For each history, see how coverage, corrected annual inputs and the conditional joint model change the simple post-protection slope. The annual-input step includes source and survey changes as well as biomass; it is not a biomass-only difference.',['series_id','Source','Taxa','Resp','old_post_years','new_post_years','shared_post_years','old_coverage_step','annual_input_step','new_coverage_step','joint_model_step']],
      'Pathway refits':['Results for review','Post-protection slopes were paired on the same years, stations and zones within each source. These conditional fits account for outcome time dependence; predictor error and cross-taxon covariance remain unresolved.',['panel','Source','candidate_MPAs','fitted_MPAs','status','slope','ci_low','ci_high','p','q']],
      'Heatwave refits':['Results for review','Updated annual ratios enter a period model with calendar time, source effects and within-history time dependence. These conditional refits do not replace the paper’s resilience analyses.',['Taxa','Resp','term','status','estimate','low','high','p','q']],
      'Open questions':['Remaining work','Questions the files do not yet answer, with the next step for each. Review the updated pathway and heatwave refits before changing the paper.',['scope','issue','finding','next_action','who_or_where']],
      'Uncertainty limits':['Remaining work','Which sources of uncertainty are included and which the available data cannot measure.',['component','what_is_included','remaining_limit']]}
    base=Path('outputs/integrated_review_20261103')
    records=lambda path: json.loads(pd.read_csv(path).to_json(orient='records'))
    visual={'comparisons':payload['Eligibility decisions'],
            'annual':records(base/'annual_values.csv'),
            'metadata':records(Path('outputs/dataset_inventory_20261101/mpa_metadata.csv')),
            'pooled':records(base/'pooled_results.csv'),
            'biomass':records(Path('outputs/biomass_policy_review_20261103/comparison_support.csv')),
            'sites':records(Path('outputs/methods_results_reference_20261104/site_decisions.csv'))}
    page=Path('code/templates/current_decisions.html').read_text()
    page=page.replace('__SITE_MAP_SCRIPT__',Path('code/templates/site_map.js').read_text())
    page=page.replace('__SITE_PAYLOAD__',Path('outputs/dataset_inventory_20261101/site_map_payload.json').read_text().replace('<','\\u003c'))
    page=page.replace('__VISUAL_SCRIPT__',Path('code/templates/decision_visuals.js').read_text())
    page=page.replace('__SITE_DIRECTORY_SCRIPT__',Path('code/templates/site_directory.js').read_text())
    page=page.replace('__VIZDATA__',json.dumps(visual).replace('</','<\\/'))
    for key,value in [('DATA',sections),('META',meta)]:
        page=page.replace('__'+key+'__',json.dumps(value).replace('</','<\\/'))
    return page
