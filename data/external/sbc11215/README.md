# SBC 112.15 Calibration Source

This is the pinned calibration source used by production harmonization from
12 October 2026. Original provider bytes are unchanged. Model 1 supplies the
recorded-month slope before temporal averaging; no pooled mean is substituted.

Rassweiler, A. A., Harrer, S., Reed, D. C., Nelson, C. J., & Miller, R. J.
(2026). *SBC LTER: REEF: Net primary production, growth and standing crop of
Macrocystis pyrifera in Southern California*, version 15. Environmental
Data Initiative. https://doi.org/10.6073/pasta/710844a10e2fe8873a4b49d97999d549

Published 9 September 2026; retrieved 12 October 2026. The package is
`knb-lter-sbc.112.15`; Table 5 is the fifth data entity, not model 5 within
that entity. Its 108 records include total foliar biomass, biomass components
and production relationships. Model 1 supplies the twelve monthly total dry
foliar standing-crop slopes, with zero intercepts.

- `metadata.xml`: exact-release EML. SHA256
  `b2369950e40c82d0672baef1aa2c0d34afb99c52c8694ae66873295713e02f32`.
- `table5.csv`: entity `6bc968141edf7b9e0a6bac923f5f78b2`, original filename
  `Regression parameters to estimate biomass and production from frond density_2026-09-09.csv`.
  Provider MD5 `50eda5d3664bf70780b3bead9674fac1`; SHA256
  `98759bc0e6825f50a82a1e077bc780be8c490991e6106996ee5f96e3ee38a0b2`.

Direct PASTA access returned HTTP 403. Both objects were recovered through
DataONE's resolve service using their exact PASTA object URLs:

- https://pasta.lternet.edu/package/metadata/eml/knb-lter-sbc/112/15
- https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/112/15/6bc968141edf7b9e0a6bac923f5f78b2

The EML releases the package under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attribute the
creators and package, retain this license information and indicate any changes.

Run `python3 -B code/audit_kelp_calibration.py` to verify bytes, identity,
quantity, units, all twelve months, the functional R conversion and the adopted
harmonized release fingerprints. See `docs/MONTHLY_KELP_REFIT_2026-10-12.md` for
the implementation and refit, and `docs/KELP_CALIBRATION_QA_2026-10-12.md` for
the historical source-verification audit before adoption.
