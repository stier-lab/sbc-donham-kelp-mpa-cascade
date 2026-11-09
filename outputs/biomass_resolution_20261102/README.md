# Biomass support: corrections and red-urchin checks

The review pipeline now withholds the Point Vicente kelp biomass comparison in 2023 because two reference-site plants lack stipe measurements. Their biomass is unknown, not zero. The whole paired year is withheld to avoid silently changing reference weights. Raw counts are retained. This produces a separate corrected review input (3,224 annual rows); the manuscript and frozen results are unchanged.

The previously checked kelp recovery slope changes from 0.243 to 0.342 log-ratio units per year; both uncertainty intervals include zero. The upstream PISCO processor now also preserves the missing-stipe flags and invalidates affected station and reference-role biomass. Its exact guard expressions were tested against the real raw and reconstructed station records. A full production rebuild remains on hold with the paper.

## Sparse sizes and borrowed sizes

We tested all 35 estimable pooled-size animal biomass recovery histories. Two additional pooled-size histories have too few post-protection years for these slope checks; two LTER lobster histories use direct lengths and are outside this pooled-size calculation. No researcher exclusion was reversed.

Shared source records receive the same random weight wherever they are used. The 500 individual-weight draws and separate 500 recorded-group draws describe how size sampling alone affects slope estimates. They are not new confidence intervals: count uncertainty, allometric uncertainty and donor representativeness are not included. Recorded group identifiers can merge multiple visits when metadata are incomplete. Existing offsets and observed years are held fixed.

A single measured animal cannot reveal its population's size variation. We therefore also show a deliberately broad stress test: halve or double singleton mean masses, with opposite changes in protected and reference pools. This is an assumption check, not a probability interval or exclusion rule. Campus Point PISCO red-urchin slopes span −0.048 to 0.040 under this test, compared with −0.005 originally. Matlahuayl red-urchin slopes span −0.392 to −0.117, compared with −0.238. These histories need caution even when a conditional bootstrap looks precise.

Keeping only donors from actual count sites retains 298 of 461 post-protection comparison-years across the 35 histories. Keeping only donors from the same program retains 300. Both restrictions permit 23 slope fits. Among those fits, the largest slope change attributable to donor choice on identical years is 0.0241 log-ratio units per year. This does not establish that omitted years or borrowed donors are harmless. The separated_changes table reports the effect of losing years separately from changing donors.

## Red urchins

All 17 red-urchin biomass histories have a dedicated role-year table and decision table. The existing equation is mass in grams = 0.00059 × diameter in millimeters^2.917. PISCO pools retain diameters of at least 25 mm; KFM and LTER use unfiltered sizes. The latter therefore represent a different measured size range. We have not removed small animals merely to make the numbers agree: any standardization must match BOTH counts and sizes.

Campus Point PISCO and LTER reuse singleton size measurements in five protected years. Those measurements do not constitute independent replication. Several island PISCO histories borrow heavily from KFM; South Point PISCO has no donors from its own count sites. The donor tables identify these explicitly instead of treating an available annual size pool as proof of local support.

Missing KFM dates remain unresolved in the combined delivery. The legacy archive contains survey dates, but linking them safely requires a verified species-code/unit crosswalk and occurrence-level matching. Matching just year, station and diameter could confuse different animals. No dates were guessed or copied into the source data.

## Recommended handling

- Use the safeguarded PISCO processor when building the next production release; use the corrected review input meanwhile.
- Retain valid density histories. Keep borrowed-size biomass provisional until donor habitat, timing and size-range compatibility are justified.
- Keep sparse-pool flags visible; do not adopt an arbitrary minimum sample cutoff or exclude by result direction.
- Carry original donor identity into joint uncertainty calculations. Complete count-and-size propagation and KFM kelp plant-sampling checks before freezing biomass uncertainty.
- For red urchins, resolve KFM date/code/unit matching and confirm that each count protocol covers the same sizes as its donor pool. The existing PISCO ≥25-mm restriction should not be imposed on other programs' size samples alone.

## Figure

Red-urchin recovery slopes on identical years: open circles use the current pooled donors; teal diamonds use only donors from actual count sites. Connecting lines show the change from donor choice. Only estimable matched-year fits appear. These are point estimates, not uncertainty intervals. Missing histories and lost years are listed in the coverage table.

## Verification

Independent R calculations reproduce 35 original slopes and uncertainty estimates, plus 46 matched-donor refits. All 487 protected manuscript/result files retain their expected hashes. The shared-donor calculation preserves 77,382 original source-row identities across 104,125 uses. The raw data were not edited. The source repository PISCO processor received the missing-stipe safeguard; its other outputs were not rebuilt. This is a review result, not a final primary analysis.
