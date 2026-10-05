// A reference identity check and an exploratory positional correspondence layer.
// Neither changes the frozen phenotype classifier or transfers an experimental allele label.
export function referenceContext(result,references){
 const ref=references.find(r=>r.id==='SNC1');
 if(!ref||ref.anchor!==637||ref.sequence.slice(636,641)!=='NLEEL')return {status:'unavailable',label:'SNC1 reference unavailable',reason:'The expected reference anchor could not be verified.'};
 const common={rule:'SNC1-context-v1',reference:'SNC1',referenceStart:637,referenceEnd:641,referenceSecondE:640,source:ref.source,experimentalEvidence:'Reported temperature-sensitive SNC1 system',conditions:ref.conditions,assay:ref.assay};
 if(result.sequence===ref.sequence)return {...common,status:'exact-reference',label:'SNC1 reference · motif verified',start:637,end:641,secondE:640,motif:'NLEEL',identity:1,coverage:1,reason:'Full query sequence exactly matches the supplied SNC1 reference. NLEEL637–641 and E640 are recognized independently of LRR-repeat numbering or the ±50-aa LRR-start cutoff. This validates reference identity and motif coordinates, not an experimental allele assignment.'};
 const comparison=result.referenceComparisons?.find(c=>c.ref.id==='SNC1'),best=comparison?.best;
 if(!best)return {...common,status:'unresolved',label:'SNC1 correspondence unresolved',reason:'No qualifying query motif window is available for reference-region alignment. A protein name alone cannot establish identity.'};
 let rp=comparison.start-1,qp=best.queryStart-1;const positions=[];
 for(let i=0;i<best.alignment.reference.length;i++){
  const a=best.alignment.reference[i],b=best.alignment.query[i];if(a!=='-')rp++;if(b!=='-')qp++;
  if(a!=='-'&&rp>=637&&rp<=641)positions.push(b==='-'?null:qp);
 }
 const contiguous=positions.length===5&&positions.every((p,i)=>p!=null&&p===positions[0]+i),coincident=contiguous&&positions[0]===best.hit.start;
 const supported=coincident&&!comparison.tied&&best.alignment.identity>=.8&&best.alignment.coverage>=.8;
 return {...common,status:supported?'aligned-candidate':'unresolved',label:supported?'Aligned SNC1-region candidate':'SNC1 correspondence unresolved',start:contiguous?positions[0]:null,end:contiguous?positions[4]:null,secondE:positions[3]??null,motif:contiguous?result.sequence.slice(positions[0]-1,positions[4]):null,identity:best.alignment.identity,coverage:best.alignment.coverage,tied:comparison.tied,reason:supported?'The ±50-aa reference motif window aligns to a qualifying query motif with ≥80% identity and paired coverage, without a tied best placement. This exploratory reference context is assessed separately from LRR-start proximity; it does not establish functional equivalence or transfer temperature sensitivity.':'The reference-window alignment is weak, tied, gapped across the motif, or does not coincide with the detected motif. Positional correspondence is not established.'};
}
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export function referenceContextHtml(r){
 const c=r.referenceContext;if(!c)return '';
 return `<section class="reference-context"><h3>${esc(c.label)}</h3><p>${esc(c.reason)}</p>${c.start?`<p><strong>Mapped reference motif: ${esc(c.motif)} ${c.start}–${c.end}; E640-equivalent: ${esc(r.sequence[c.secondE-1])}${c.secondE}.</strong> LRR-region start: ${r.domainEvidence.lrrRegionStart??'unknown'}; offset: ${r.domainEvidence.lrrRegionStart==null?'unknown':c.start-r.domainEvidence.lrrRegionStart} aa. Repeat containment and proximity remain reported separately.</p>`:''}${c.status==='exact-reference'?`<p><strong>${esc(c.experimentalEvidence)}</strong> · ${esc(c.conditions)} ${c.source?`<a href="${esc(c.source)}" target="_blank" rel="noopener">Published source</a>`:''}</p><p class="tiny">Reference-context validation: passed. Frozen sequence-only prediction: ${esc(r.classification.label)}. The experimental allele/assay context remains separate; no phenotype-ledger verification flag is changed.</p>`:'<p class="tiny">Exploratory correspondence only. The 80% window criteria are inspection thresholds, not calibrated phenotype probabilities. No phenotype is transferred from SNC1.</p>'}</section>`;
}
