import {parseDomains} from './domains.mjs';
const BASE='https://www.ebi.ac.uk/Tools/services/rest/iprscan5';
export function domainKind(d){const s=`${d.name} ${d.accession} ${d.interproName??''}`;return /NB.ARC|PF00931|IPR002182/i.test(s)?'NB-ARC':/leucine.rich.repeat|\bLRR(?:_|\b)|SSF52058/i.test(s)?'LRR':/kinase/i.test(s)?'Kinase':/\bTIR(?:_|\b)|Toll.*interleukin/i.test(s)?'TIR':/winged.helix|SSF46785/i.test(s)?'WHD candidate':'Other';}
export function deriveRegionAnnotation(annotation){
 const domains=(annotation.domains??[]).map(d=>({...d,kind:domainKind(d)}));
 const result={...annotation,domains};const lrr=domains.filter(d=>d.kind==='LRR');
 // Choose one database's span. Retain all individual calls and alternatives for inspection.
 const preferred=lrr.filter(d=>/Pfam/i.test(d.source));const spans=preferred.length?preferred:lrr;
 if(spans.length){result.lrrRegionStart=Math.min(...spans.map(d=>d.start));result.lrrRegionEnd=Math.max(...spans.map(d=>d.end));result.lrrRegionSource=`${[...new Set(spans.map(d=>d.source))].join('; ')}; envelope of LRR signature spans (gaps may be unannotated)`;}
 const nb=domains.filter(d=>d.kind==='NB-ARC');
 const whd=domains.filter(d=>d.kind==='WHD candidate'&&nb.some(n=>d.start<=n.end&&d.end>=n.start)&&d.end<(result.lrrRegionStart??result.firstLrrStart??Infinity));
 if(whd.length===1&&!result.whdStart){result.whdStart=whd[0].start;result.whdEnd=whd[0].end;}
 return result;
}
async function request(path,options={},fetcher=fetch){let r;try{r=await fetcher(`${BASE}/${path}`,{...options,signal:options.signal??AbortSignal.timeout(60000)});}catch(e){throw new Error(`InterProScan connection interrupted. ${e.message}`);}const text=await r.text();if(!r.ok)throw new Error(`InterProScan HTTP ${r.status}: ${text.replace(/<[^>]*>/g,' ').slice(0,180)}`);return text;}
export async function submitDomainJob(record,email,fetcher=fetch){
 if(!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email))throw new Error('Enter a valid job email for InterProScan.');
 if(!/^[ACDEFGHIKLMNPQRSTVWYXBZJUO]+$/.test(record.sequence))throw new Error('Invalid protein sequence.');
 const body=new URLSearchParams({email,stype:'p',sequence:`>nlr_query\n${record.sequence}`,title:'NLR-Temper domain annotation',goterms:'false',pathways:'false'});
 const id=(await request('run',{method:'POST',body},fetcher)).trim();if(!/^iprscan5-[a-zA-Z0-9_-]+$/.test(id))throw new Error('InterProScan did not return a valid job identifier.');return id;
}
export async function domainJobStatus(id,fetcher=fetch){if(!/^iprscan5-[a-zA-Z0-9_-]+$/.test(id))throw new Error('Invalid job ID.');return (await request(`status/${id}`,{},fetcher)).trim();}
export async function collectDomainJob(id,record,fetcher=fetch){
 if(!/^iprscan5-[a-zA-Z0-9_-]+$/.test(id))throw new Error('Invalid job ID.');const text=await request(`result/${id}/tsv`,{},fetcher);
 const source=`InterProScan / EMBL-EBI job ${id}`;
 if(!text.trim())return {source,domains:[],domainJob:{id,status:'complete',noHits:true}};
 const a=parseDomains(text,'interpro',[{id:'nlr_query',sequence:record.sequence}]).get('nlr_query');
 if(!a)throw new Error('Returned annotations did not match the submitted query.');
 return deriveRegionAnnotation({...a,source,domainJob:{id,status:'complete',noHits:false}});
}
