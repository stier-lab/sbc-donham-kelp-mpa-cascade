"""Admission guard for review biomass draws; never changes density or raw data."""
EXCLUDED_KFM = {f'a-k-{n:02}' for n in (21,5,7,26,27,28,29,30,36,37,12,13,31,35)}
URCHINS = {'M. franciscanus','S. purpuratus'}
def source_issue(site, taxon, source_table='ALL_sizefreq_2024.csv'):
    if site in EXCLUDED_KFM:
        return 'Explicit researcher exclusion'
    if source_table == 'ALL_sizefreq_2024.csv' and taxon in URCHINS and site in {'a-k-01','a-k-14'}:
        return 'Known composite or wrong-site size label; original permitted records required'
    return ''
def require_allowed(site, taxon, source_table='ALL_sizefreq_2024.csv'):
    issue=source_issue(site,taxon,source_table)
    if issue: raise ValueError(issue+': '+site)
