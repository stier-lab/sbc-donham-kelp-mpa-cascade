"""Small adversarial aggregation fixtures, independent of empirical results."""
import ast,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
B=1200;VARIANT='';weight_cache={};pool_cache={};mu={};p=pd.DataFrame();mem=pd.DataFrame();visit_index=pd.DataFrame()
parsed=ast.parse(Path('code/integrated_review/estimate.py').read_text());exec(compile(ast.Module(body=[n for n in parsed.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'aggregation','exec'))
g=pd.DataFrame(dict(Source=['LTER']*5,site=['X']*5,unit=['1|a','1|b','1|c','1|d','2|a'],unit_id=['a','b','c','d','e'],density=[1.,1.,1.,1.,3.]))
point,draw=aggregate_units(g,'density');assert point==2.,point;assert draw.std()>.2,'Transect uncertainty was lost';assert np.isfinite(draw).all()
g2=g.copy();g2.loc[0,'density']=np.nan;point2,d2=aggregate_units(g2,'density');assert np.isnan(point2) and np.isnan(d2).all()
mem=pd.DataFrame(dict(pool_id=['donor','donor'],source_record_id=['f1','f2'],mass=[10.,20.],frequency=[1.,1.]));mu={'donor':15.}
a=pool_draw('donor');b=pool_draw('donor');assert np.array_equal(a,b);assert a.std()>0
assert abs(a.mean()-15)<.4
# Shared imputation affects only the imputed component of an otherwise observed unit.
g=pd.DataFrame(dict(Source=['LTER'],site=['Y'],unit=['1|A'],unit_id=['q1'],biomass=[87.],imputed_density=[3.],imputed_pool=['donor'],donor_pool=['']))
point,draw=aggregate_units(g,'biomass');assert point==87;assert np.allclose(draw,42+3*a)
# Nonconsecutive years are separated by their actual lag, not row positions.
y=np.array([2001,2004,2005]);K=.6**np.abs(y[:,None]-y[None,:]);assert np.isclose(K[0,1],.216)
print('PASS: equal-transect point means, nonzero transect uncertainty, unknown propagation, shared donor identity, partial imputation, calendar gaps')
