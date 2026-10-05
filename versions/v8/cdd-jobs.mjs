import {deriveRegionAnnotation} from './domain-jobs.mjs';
export function parseCddStatus(text){const status=Number(text.match(/^#status\s+(\d+)/m)?.[1]??-1),id=text.match(/^#cdsid\s+(\S+)/m)?.[1];if(![0,3].includes(status))throw new Error('NCBI CDD: '+({1:'invalid or expired job ID',2:'no effective input',4:'queue service error',5:'results expired or unavailable'}[status]??'unrecognized response'));return {id,status:status===0?'FINISHED':'RUNNING',text};}
export async function cddRequest(input,fetcher=fetch){const r=await fetcher('/api/cdd',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(input),signal:AbortSignal.timeout(60000)});const out=await r.json();if(!r.ok||out.error)throw new Error(out.error||'CDD request failed');return parseCddStatus(out.text);}
export function collectCddResult(text,record,id){
 const status=parseCddStatus(text);if(status.status!=='FINISHED')throw new Error('CDD job is still running.');const domains=[];
 for(const line of text.split(/\r?\n/)){if(!line.startsWith('Q#'))continue;const c=line.split('\t');if(c[0].startsWith('Q#N'))continue;if(!/^Q#1\s*-\s*>?nlr_query\b/.test(c[0]))throw new Error('CDD results do not match this submitted query.');
  const start=Number(c[3]),end=Number(c[4]),e=Number(c[5]);if(!Number.isInteger(start)||!Number.isInteger(end)||start<1||end<start||end>record.sequence.length||!Number.isFinite(e))throw new Error('Invalid CDD domain coordinates.');
  if(domains.some(d=>d.start===start&&d.end===end&&d.accession===c[7]&&d.hitType===c[1]))continue;
  domains.push({start,end,name:c[8],accession:c[7],source:'NCBI CDD / '+c[1],engine:'cdd',score:c[5],bitScore:Number(c[6]),pssm:c[2],hitType:c[1]});
 }
 return deriveRegionAnnotation({source:'NCBI CDD job '+id,domains,domainJob:{id,engine:'cdd',status:'complete',noHits:domains.length===0}});
}
