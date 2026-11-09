# Naples and satellite source reconciliation

**19 October update:** Jenn Caselle confirmed that KFM and PISCO always collect
separate observations, including at shared sites or nearby dates. The candidate
ledger below is a spatial/temporal proximity check, not evidence of record reuse.
Separate collection does not remove ecological dependence or shared-reference
covariance. The original report code exists but is not readily runnable; its
exact annual inputs remain unavailable. See the current sheephead YOY audit for
the supplied lookup's size threshold and the retained-record sensitivity.

## What is resolved

The current Naples PISCO references match the published monitoring report.
The remaining reference question is whether contributors intended a later
replacement, not whether the current mapping occurs in that report.
Satellite units and many annual values are verified, but the delivered
missingness and extraction version are not fully reconstructed. Nearby
same-day surveys can now be identified for contributor review; none has been
classified as a duplicate.

No primary data, model, figure or inclusion rule was changed. The essential
supplement's S1.3 now states the verified Naples mapping. Detailed diagnostics
stay in the repository rather than extending the supplement.

## Naples and the technical report

The recovered source is *Monitoring and Evaluation of Kelp Forest Ecosystems
in the MLPA Marine Protected Area Network*, dated 30 December 2021, in the
[university-hosted v2 PDF](https://kelpforest.msi.ucsb.edu/assets/files/1658/Kelp%20Forest%20Technical%20Report%20Narrative_v2.pdf).
Appendix Table 1, PDF page 116 (appendix printed page 5), lists
ARROYO_QUEMADO_E and ARROYO_QUEMADO_W as references for Naples/Campus Point.
Both rows were visually checked. This agrees with our present Naples PISCO
mapping; it does not establish that every later internal report used it.

The report's focal-fish response is biomass; focal invertebrates and algae use
density. It fits log inside/reference ratios against calendar year, whereas
our primary results combine selected mean-level and fitted-change contrasts.
The report excludes fish canopy transects and YOY fish (PDF pp. 40-41,
printed pp. 38-39). These are documented differences in analysis definitions,
not proof that any one difference caused the reported discrepancy.

Two Naples coefficients were transcribed and visually verified: sheephead
biomass, Table 19, PDF p. 167: slope 0.008532, SE 0.027163, P = 0.7599;
giant-kelp density, Table 20, PDF p. 177: slope 0.024685, SE 0.016166,
P = 0.1611. `report_naples_verified.csv` records the source locations.
Applying calendar-year regressions to current PISCO log ratios through 2020
gives P = 0.8583 (sheephead, 12 years) and P = 0.2844 (kelp, 11 years).
Neither is an exact reproduction. A change of log base rescales coefficients
but cannot explain different P values. The report's log base, zero handling,
exact annual input version and original analysis code remain unverified.

**Still ask Jenn/Emily:** Did any later analysis replace these Arroyo Quemado
references? If so, which sites, years and file version? For exact reproduction,
obtain the annual input and script behind the specific disputed figure. The
public data release associated with later work is not automatically the
2021 report's input version.

## Satellite delivery

The report also describes survey-site polygons 200 m alongshore, bounded by
the 5 m contour and an offshore transect waypoint or 20 m contour (PDF p. 29,
printed p. 27). Thus site-footprint products existed in that project. It does
not establish that our delivered whole-MPA CSV contains those products.
Its correlation-selected reference recipe, using 1984-2020 histories and
candidate regions within approximately 50 km, must not be substituted for
the later 2022 pixel workbook's selection recipe.

The diagnostic rebuild uses the recovered 2022 workbook and linked NetCDF,
summing finite wet-biomass pixel values by quarter before annual maxima.
Among 60 delivered rows, 38 have an exact candidate match on all paired finite
years; only six match the entire annual vector, including missingness. The
other 32 rows omit 46 row-years for which the recovered pixel file yields a
finite annual maximum. Most are in 1984, but omissions also occur in 1988-1993
and 2012. Some omitted years have four finite quarterly sums. This rules out
an across-the-board incomplete-year rule as the sole explanation. Filling
all-missing reconstructed years with zero closes none of these mismatches.
Neither calculation justifies filling the delivered gaps.

For 31 nonempty reference sets with a mapped MPA, 21 reference centroids are
outside 1-3 km of the mean MPA-pixel coordinate. This tests centroid-to-centroid
distance, not merely the outermost reference pixel. No reference set shares
pixels with its own workbook MPA set. This does not verify exclusion from
every other MPA. Exact spherical-polygon centroids, source mask versions and
any exceptions to the emailed recipe still require the extraction script.

**Still ask Tom/Emily:** Which script and masks produced the delivered CSV;
what rule removed these specific row-years; which reference-selection version
applies; and where are the common-survey-footprint polygons and extraction?
The emitted `satellite_missing_years.csv` makes the missingness question
specific. Do not silently swap reference replicate numbers: several workbook
groups match the opposite delivered replicate.

## Survey overlap

The MBON files contain PISCO records as well as KFM, but the KFM processor
explicitly selects `data_source == "kfm"` for benthos and fish; LTER is split
separately. Primary KFM sheephead uses `visualfish`, not roving counts.
Consequently, the mere presence of PISCO rows in the MBON download is not
evidence that our pipeline reused them as KFM.

The source audit streams all supplied MBON and PISCO fish/swath inputs,
records source/method counts and SHA256 hashes, and compares valid survey
dates between KFM and PISCO site centres within 1 km. The ledger contains
175 site-by-survey-kind coordinate pairs, of which 17 have at least one
same-day visit. It deliberately covers more than the retained primary cohort.
Repeated site-to-MPA mappings are collapsed only for coordinate lookup;
observations are not deduplicated. Some nominal sites have multiple
coordinates, so the ledger is not a count of unique reefs.

Same-day proximity does not establish shared observations: methods, areas,
transect identities and source record identifiers must be reconciled first.
Nor do separate dates establish ecological independence. Programme-separated
matched-year refits answer a different question and do not resolve those
identities. The current supplement retains that limitation.

**Still ask Jenn/programme leads:** Are the listed nearby same-day surveys
independent collections, shared field visits with distinct measurements, or
reused observations? Use the candidate dates rather than a general question
about programme overlap. No rows were removed on proximity alone.

## Rebuild and scope

Outputs are under `data/external/source_reconciliation_2026-10-17/`.
Original report/workbook/NetCDF files remain in ignored local storage.

```sh
python code/tools/audit_survey_sources.py /path/to/data-repository --report /path/to/report-v2.pdf
Rscript code/tools/reconcile_satellite_delivery.R /path/to/CAkelpCanopyEnv_2021.nc
Rscript code/tools/compare_naples_report.R
```

These are read-only source diagnostics, not pipeline steps that alter primary
measurements. Scope agreement about Sofia's work still needs Jenn and the
affected collaborators; file comparisons cannot provide that agreement.
No new clarification email was sent in this follow-up.
