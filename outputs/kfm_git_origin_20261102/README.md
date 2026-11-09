**November 2 follow-up correction:** inspection of Emily’s saved R workspace and original site table shows that KFM sizes were not universally unmatched before the June backfill. Wyckoff already has a matching Harris Point reference key. The June change cannot be described as the universal activation of these errors. See the [original-folder investigation](../kfm_pipeline_guard_20261102/README.md) and saved-workspace membership table for the narrower evidence.

# How the KFM site-label errors reached biomass estimates

## Finding

There are two separate failures: the combined input already contains misleading physical-site labels; later processing trusts those labels when assigning sizes to an MPA. Git identifies the downstream code change, but not the original creation of the bad labels. Do not attribute that upstream recoding to a person or model without the missing assembly code or dated source revisions.

## Verified history

- The data repository's initial February 27, 2026 commit (`fd3045e2f98d334ca1f5d2a6c1aad21f14f67d9f`) reads `ALL_sizefreq_2024.csv` and excludes delivered site IDs, including a-k-21 and a-k-37. It does not build the combined size file. The retained original V10 script likewise reads it. Neither checked script creates the two bad labels.
- June 17 commit [93a02d3](https://github.com/stier-lab/sbc-donham-kelp-mpa-data/commit/93a02d347136e914aadb3b3b6fb1da7453b74667) adds MPA/status metadata to previously unmatched KFM/LTER size rows using their site IDs. That enables these size rows to enter the biomass pools. The code takes the delivered IDs as correct; it does not check original physical-site identity. This change is recorded under Adrian's Git identity with a Claude coauthor trailer. That identifies the recorded implementation, not who originally assembled the erroneous input.
- The current source is a symlink to the shared Drive file. The preserved original-project copy is byte-identical. A committed October 7 snapshot in `adrianstier/sbc-sswd-kelp` (`a70c06397a5c5972cc9afe15debd24dae6e91a9b`) has exactly the same KFM site/year/species/size frequencies. This proves the input problem existed by October 7. The raw combined file is absent from the data/analysis repositories' tracked histories, so those histories do not prove its exact state on June 17 or when it was first mislabeled.
- November 2 reconciliation establishes the physical identity by matching every exact size and frequency against Scott's original export. Wyckoff's delivered label combines Wyckoff and Miracle Mile; the delivered Southeast Sea Lion label contains Southeast Reef. Excluding IDs21/37 cannot remove measurements already labeled01/14.

The June code change explains the route into biomass, not the origin of the bad labels. Whether the upstream assembly used a wrong lookup, manual recoding, or another operation remains unknown. No code evidence supports choosing among those explanations.

## Broader source-identity screen

Compared all KFM red/purple urchin distributions in the combined file for years represented in Scott's export (2012–2019 and 2021–2022), not only the six affected response histories. Exact frequency comparisons identify a unique match for all558site/species/year groups across28delivered site labels:

- 518 groups across26labels match their intended physical sites. Gull Island/Gull Island South and Black Seabass/Black Sea Bass are naming differences, with matching distributions.
- 20groups under a-k-01 exactly match Wyckoff plus Miracle Mile.
- 20groups under a-k-14 exactly match Southeast Reef, not Southeast Sea Lion.
- No unmatched group and no additional wrong-site distribution found within this scope. All two-site combinations were tested when no single-site distribution matched.

This is not a clean bill of health for every dataset. Scott's export contains only red/purple urchins, and this broad screen covers its available years. Other taxa and earlier years require their own original-source checks. Older source rows remain a broader audit target beyond the already verified historical corrections.

## Other related weakness

The MBON site table has three legitimate reference sites (Admiral's Reef, East Fish Camp, Lighthouse) assigned to three Anacapa comparison histories. A plain first-match lookup cannot represent all those memberships; it can silently retain only one. This is a comparison-membership issue, not proof of a wrong physical site. The earlier November1 [Anacapa audit](../upstream_sites_20261101/README.md) identified and repaired SMCA membership in isolated review data. The legacy production backfill still uses `match()` and must not be treated as a fully repaired primary pipeline. Separate physical identity from comparison membership and use an explicit site–comparison key.

## Prevention and remaining work

- Preserve original physical site ID/name on each size measurement before harmonization. Blank original KFM names in the combined file remove this cross-check.
- Apply researcher exclusions to physical IDs before pooling; retain excluded-record counts and source lineage.
- Require a one-to-one source-to-physical-site mapping. Explicitly declare many-comparison reuse; never resolve it by table order.
- Validate source frequency totals, dates, units and species before each biomass rebuild. The exact-distribution audit here is a reproducible check, not just a visual label review.
- Promote the verified replacements through one authoritative ingestion path before rebuilding primary results; currently they are review-only corrections.
- Locate the script/export that originally assembled `ALL_sizefreq_2024.csv` and its site lookup. The checked Git histories, original-project R files/R history and GitHub code search did not locate it. This is the specific missing evidence for the initial mistake; no new email was sent.

## Evidence and reproduction

Run `python3 code/tools/audit_kfm_site_identity_history.py`. Outputs: [all558matches](all_kfm_site_year_matches.csv), [site summary](site_identity_summary.csv), [source-copy hashes](source_copy_comparison.csv), [dated snapshot](dated_snapshot.json), [GitHub commit metadata](github_commit.json), [exact code patch](metadata_activation.patch), [multiple comparison memberships](multiple_comparison_memberships.csv). Source measurements remain in their existing private/raw locations. Analysis and data remotes were fetched without modifying working files. No manuscript, primary numerical result, source data or coauthor communication changed.
