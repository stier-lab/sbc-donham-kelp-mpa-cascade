const fs=require('fs'),vm=require('vm'),assert=require('assert');
const nodes={};const ctx=vm.createContext({document:{getElementById(id){return nodes[id]??={innerHTML:'',textContent:'',value:'',checked:false,scrollIntoView(){}}}}});
const page=fs.readFileSync('outputs/dataset_inventory_20261101/index.html','utf8');
const scripts=[...page.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(x=>x[1]);
vm.runInContext(scripts.find(x=>x.includes('const data=')),ctx);
const run=s=>vm.runInContext(s,ctx);
assert.equal(run('cmp.length'),206);assert.equal(run('cmp.filter(c=>c.exclude).length'),31);assert.equal(run('models.length'),175);
assert(run('models.every(m=>String(m.years||"").split(";").filter(Boolean).every(y=>yearState(m.series_id,Number(y))==="supported"))'));
assert(run('viz.annual.filter(r=>r.response_supported).every(r=>String(modelById.get(r.series_id).years).split(";").includes(String(r.year)))'));
assert(run('roster.filter(r=>!r.survey_eligible).every(r=>yearState(r.series_id,r.year)==="coverage")'));
assert(run('cmp.filter(c=>c.exclude).every(c=>!modelById.has(c.series_id))'));
assert(nodes.decisionFlow.innerHTML.includes('62 before–after'));assert(nodes.decisionFlow.innerHTML.includes('100 later trend only'));
assert((nodes.pooledSnapshot.innerHTML.match(/<tr>/g)||[]).length===10);
assert(nodes.pooledSnapshot.innerHTML.includes('1449 of 1990'));
assert(nodes.pooledSnapshot.innerHTML.includes('20 PISCO comparison-years'));
for(const site of run('[...new Set(cmp.map(c=>c.MPA))]')){
 run('selectSite('+JSON.stringify(site)+')');assert(nodes.timeline.innerHTML.includes('Program · species'));assert(nodes.comparisonDetails.innerHTML.includes(site));
}
run('selectSite("Harris Point SMR")');assert(nodes.stationMap.innerHTML.includes('<path'));assert(nodes.stationMap.innerHTML.includes('<rect'));assert(nodes.mapList.innerHTML.includes('Hare Rock'));assert(nodes.siteMapSection.hidden===false);
run('selectSite("Abalone Cove SMCA")');assert(nodes.siteMapSection.hidden===true);
run('selectSite("Naples SMCA")');run('pickComparison(selectedIndex,2012)');assert(nodes.comparisonDetails.innerHTML.includes('2012'));assert(nodes.timeline.innerHTML.includes('designation year'));
run('tableLink("Survey coverage",cmp[selectedIndex].series_id)');assert(nodes.count.textContent!=='0 rows');
fs.writeFileSync('outputs/dataset_inventory_20261101/visual_checks.json',JSON.stringify({candidate_histories:206,excluded:31,retained:175,before_after:62,post_only:100,insufficient_post:13,supported_years_agree_with_models:true,all_sites_render:true,coverage_failures_distinct:true,year_selection:true,table_links:true},null,2)+'\n');
console.log('Visual state agrees with current annual values, model eligibility and exclusions; all sites, year selection and table links passed.');
