# Seasonal coverage of the retained comparisons

Seasonal growth, recruitment and movement could produce protected–reference
differences when the reefs are surveyed in different months. Annual matching
alone does not remove this possibility.

Module 54 reconstructs recorded count-survey months for every retained diver
history and year, using the source-specific station, taxon and survey filters.
It retains surveyed zeros, separates density and biomass histories, and audits
separate kelp and animal size samples. Unique site-months prevent species-rich
or heavily replicated visits from dominating the coverage counts. These are
survey-coverage bounds: a visit does not prove every biomass record was valid.

Across retained response-years, the two sides share at least one survey month
in 1,025/1,260 PISCO, 474/800 KFM and 196/198 LTER comparisons. The remaining
235, 326 and 2 response-years, respectively, have no common survey month.
These totals repeat a survey when it contributes to multiple taxa/responses;
they are not independent sample sizes. Common months do not imply identical
dates, equal survey effort, or that an annual mean uses only the common months.

PISCO benthic surveys use May–October records; retained fish survey coverage
extends outside that window. Retained KFM count surveys span May–October,
whereas retained LTER surveys span July–October. Thus handbook windows do not
fully describe the actual inputs. KFM-derived urchin size pools lack months in
the supplied combined size file; their seasonal agreement with counts cannot
be established from that file. Known size dates remain separate from counts
in the audit, including KFM kelp's dated plant samples.

The audit establishes a seasonal limitation, not its direction or magnitude.
No annual abundance estimate or selected model is changed by this audit.
A same-season abundance sensitivity is warranted, particularly for KFM and
PISCO; it must reconstruct responses from common-month surveys, retain the same
sites and years where possible, and compare the original and restricted methods
on that common coverage. Dropping years with disjoint month sets alone would
confound seasonal matching with a change in temporal coverage. Biomass checks
also need dated, compatible size samples. Neither month-specific kelp mass
conversion nor annual covariance corrections removes survey-timing differences.

Run `Rscript code/R/54_seasonal_coverage.R` then
`python3 code/tests/check_seasonal_coverage.py`. Override `DONHAM_DATA_REPO` for
the source-data workspace. Tables `table_s_seasonal_*` record site-months,
paired-year overlap, period summaries, separate size dates and input hashes.
The main text states the timing limitation; SI S1.1 gives the program totals
and points to these detailed repository tables.

## October 30 user scope decision

Adrian explicitly deferred seasonal reanalysis and sensitivities to another time. Retain the existing brief main-text disclosure and completed audit/SI support. Seasonal refitting is not an active task or a blocker for this revision; Li’s seasonal disclosure request is addressed. No primary analysis or manuscript files changed for this scope decision.
