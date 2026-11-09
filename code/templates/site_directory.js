const siteRows=viz.sites.map((row,index)=>({...row,index}));
const siteComparisonIndex=new Map(cmp.map((row,index)=>[row.series_id,index]));
const directoryPageSize=24;
let siteDirectoryArea=selectedSite,siteDirectoryProgram='',siteDirectoryQuery='',siteDirectoryPage=0,siteDirectoryIndex=null,siteDirectoryReady=false;

function siteTitle(row){return row.site_name||row.site_id||'Unnamed source record'}
function siteRole(role){return role==='mpa'?'Protected':role==='reference'?'Reference':role==='not_assigned'?'No current role':String(role||'Role unknown').replaceAll('mpa','protected')}
function siteStatus(status){return {
  included_for_some_responses:'Used for at least one species or measurement',
  excluded_from_retained_responses:'Assigned, but not used for the linked measurements',
  not_in_current_comparison_roster:'Not assigned to a current comparison'
}[status]||'Status needs checking'}
function siteSeriesIds(value){return String(value||'').split(';').filter(Boolean)}
function siteMatches(row){
  if(siteDirectoryArea&&row.MPA!==siteDirectoryArea)return false;
  if(siteDirectoryProgram&&row.program!==siteDirectoryProgram)return false;
  if(!siteDirectoryQuery)return true;
  return [row.site_name,row.site_id,row.MPA,row.program,row.role,row.included_responses,row.excluded_responses].some(value=>String(value||'').toLowerCase().includes(siteDirectoryQuery));
}
function siteSeriesLinks(ids){
  return ids.map(id=>{
    const i=siteComparisonIndex.get(id),comparison=cmp[i];
    if(i===undefined||!comparison)return '';
    return `<button type="button" onclick="openSiteComparison(${i})">${esc(seriesLabel(comparison))}</button>`;
  }).filter(Boolean).join(' · ')||'None recorded';
}
function renderSiteRecord(){
  const row=siteRows[siteDirectoryIndex],card=document.getElementById('directoryCard');
  if(!row){card.innerHTML='<p>No site matches these filters.</p>';return}
  const first=Number(row.first_raw_year),last=Number(row.last_raw_year);
  const span=Number.isFinite(first)&&Number.isFinite(last)&&first>0&&last>0?(first===last?String(first):`${first}–${last}`):'Years not recorded';
  const min=Number(row.depth_min_m),max=Number(row.depth_max_m);
  const depth=row.depth_min_m!==null&&row.depth_max_m!==null&&Number.isFinite(min)&&Number.isFinite(max)?`${min}–${max} m across source records`:'No comparable depth range recorded';
  const included=siteSeriesIds(row.included_series_ids),excluded=siteSeriesIds(row.excluded_series_ids);
  const reason=row.decision_reasons?`<p><b>Current assignment note:</b> ${esc(row.decision_reasons)}</p>`:'';
  const size=row.size_status?`<p><b>Size data:</b> ${esc(row.size_status)}</p>`:'';
  card.innerHTML=`<h3>${esc(siteTitle(row))}</h3>${row.site_name?'':'<p class="caption">No separate site name is recorded; the site code is shown above.</p>'}<p><span class="site-status">${esc(siteStatus(row.roster_status))}</span></p><p><b>Site code:</b> ${esc(row.site_id||'Not recorded')} · <b>Program:</b> ${esc(row.program)} · <b>Area:</b> ${esc(row.MPA||'Not recorded')}</p><p><b>Role:</b> ${esc(siteRole(row.role))} · <b>Years in source:</b> ${esc(span)}</p><p><b>Recorded depth:</b> ${esc(depth)}</p>${reason}${size}<p><b>Used for:</b> ${siteSeriesLinks(included)}</p><p><b>Assigned but left out for:</b> ${siteSeriesLinks(excluded)}</p>${!included.length&&!excluded.length?'<p>This record has no assignment to a current comparison. That alone does not mean a survey was absent or that the site was rejected for ecological reasons.</p>':''}<p class="caption">Source years and depth describe this site record; usable years and shared survey zones are checked separately for each comparison.</p>`;
}
function renderSiteDirectory(){
  const matches=siteRows.filter(siteMatches).sort((a,b)=>String(a.program||'').localeCompare(String(b.program||''))||siteTitle(a).localeCompare(siteTitle(b))||String(a.site_id||'').localeCompare(String(b.site_id||'')));
  const pages=Math.max(1,Math.ceil(matches.length/directoryPageSize));siteDirectoryPage=Math.max(0,Math.min(siteDirectoryPage,pages-1));
  const shown=matches.slice(siteDirectoryPage*directoryPageSize,(siteDirectoryPage+1)*directoryPageSize);
  if(!shown.some(row=>row.index===siteDirectoryIndex))siteDirectoryIndex=shown[0]?.index??null;
  document.getElementById('directoryCount').textContent=`${matches.length} site or group ${matches.length===1?'record matches':'records match'}. ${matches.filter(row=>row.roster_status==='included_for_some_responses').length} ${matches.filter(row=>row.roster_status==='included_for_some_responses').length===1?'is':'are'} used for at least one measurement.`;
  document.getElementById('directoryList').innerHTML=shown.map(row=>`<button type="button" aria-current="${row.index===siteDirectoryIndex}" onclick="pickSiteRecord(${row.index})">${esc(siteTitle(row))}<small>${esc(row.program)} · ${esc(siteRole(row.role))} · ${esc(row.site_id||'No site code')}</small></button>`).join('')||'<p class="caption">No matching sites.</p>';
  document.getElementById('directoryPage').textContent=`Page ${siteDirectoryPage+1} of ${pages}`;
  document.getElementById('directoryPrev').disabled=siteDirectoryPage===0;
  document.getElementById('directoryNext').disabled=siteDirectoryPage>=pages-1;
  renderSiteRecord();
}
function pickSiteRecord(index){siteDirectoryIndex=index;const area=siteRows[index].MPA;if(cmp.some(c=>c.MPA===area))selectSite(area,false);renderSiteDirectory()}
function moveSiteDirectory(delta){siteDirectoryPage+=delta;siteDirectoryIndex=null;renderSiteDirectory()}
function openSiteComparison(index){const comparison=cmp[index];if(!comparison)return;selectSite(comparison.MPA);pickComparison(index);document.getElementById('comparisonDetails').scrollIntoView({behavior:'smooth'})}
function syncSiteDirectory(area){siteDirectoryArea=area;siteDirectoryProgram='';siteDirectoryQuery='';siteDirectoryPage=0;siteDirectoryIndex=null;document.getElementById('directoryArea').value=area;document.getElementById('directoryProgram').value='';document.getElementById('directorySearch').value='';renderSiteDirectory()}
function initSiteDirectory(){
  const area=document.getElementById('directoryArea'),program=document.getElementById('directoryProgram'),search=document.getElementById('directorySearch');
  area.innerHTML='<option value="">All areas</option>'+[...new Set(siteRows.map(row=>row.MPA))].sort().map(value=>`<option value="${esc(value)}">${esc(value==='[blank]'?'Area not recorded':value)}</option>`).join('');
  program.innerHTML='<option value="">All programs</option>'+[...new Set(siteRows.map(row=>row.program))].sort().map(value=>`<option value="${esc(value)}">${esc(value)}</option>`).join('');
  area.value=siteDirectoryArea;
  area.onchange=()=>{siteDirectoryArea=area.value;siteDirectoryPage=0;siteDirectoryIndex=null;if(area.value&&cmp.some(c=>c.MPA===area.value))selectSite(area.value,false);renderSiteDirectory()};
  program.onchange=()=>{siteDirectoryProgram=program.value;siteDirectoryPage=0;siteDirectoryIndex=null;renderSiteDirectory()};
  search.oninput=()=>{siteDirectoryQuery=search.value.trim().toLowerCase();if(siteDirectoryQuery){siteDirectoryArea='';siteDirectoryProgram='';area.value='';program.value=''}siteDirectoryPage=0;siteDirectoryIndex=null;renderSiteDirectory()};
  document.getElementById('directoryPrev').onclick=()=>moveSiteDirectory(-1);
  document.getElementById('directoryNext').onclick=()=>moveSiteDirectory(1);
  siteDirectoryReady=true;renderSiteDirectory();
}
