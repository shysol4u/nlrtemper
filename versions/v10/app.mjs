import {jointReferenceComparison,jointReferenceHtml} from './joint-reference.mjs';
import {resolveStructureSelection,authorCoordinate} from './structure-selection.mjs';
import {shortReferenceComparisons,shortReferenceHtml} from './short-reference.mjs';
import {cddRequest,collectCddResult} from './cdd-jobs.mjs';
import {deriveRegionAnnotation,submitDomainJob,domainJobStatus,collectDomainJob} from './domain-jobs.mjs';
import {residueContext} from './residue-context.mjs';
import {parseDomains} from './domains.mjs';
import {referenceComparisons} from './reference.mjs';
import {architecture,numberedSequence,referenceHtml} from './visuals.mjs';
import {analyzeRecord,contextFor,parseFasta,integratedCall} from './core.mjs';
import {EXAMPLES} from './examples.mjs';
import {createStructureViewer,fetchAlphaFoldPdb,inferUniProtAccession} from './structure.mjs';
import {resultsTable,delimited} from './export.mjs';
const $=s=>document.querySelector(s),esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'}[c]));
const pct=v=>v==null?'—':`${Math.round(v*100)}%`,loc=m=>m?`${m.sequence} · ${m.start}–${m.end}`:'None';
let references=[],selectedHits=new Map(),domainJobs=new Map(),domainBusy=false,stopDomainJobs=false;
let structureChoices=new Map();
let config,ledger={records:[]},report=null,results=[],annotations=new Map(),viewer,chains=[],payload=null,structureSequence=null;
const badge=label=>`<span class="badge ${label==='TS-like'?'ts':label==='TT-like'?'tt':'uncertain'}">${esc(label)}</span>`;
const safeLink=url=>/^https:\/\//.test(url)?`<a href="${esc(url)}" target="_blank" rel="noopener">Source</a>`:'Source pending';
const phenotypeRecords=r=>ledger.records.filter(p=>p.sequence===r.sequence);
function combine(r){r.jointComparison=jointReferenceComparison(r,references);r.shortComparisons=shortReferenceComparisons(r,references);r.referenceComparisons=referenceComparisons(r,references);r.integrated=integratedCall(r,phenotypeRecords(r),r.structure);return r;}
function sequenceHtml(sequence,hits){return numberedSequence(sequence,hits);}
function schematic(r){return architecture(r,results.indexOf(r));}
function hitTable(r){return `<div class="table-wrap"><table><thead><tr><th>Sequence</th><th>Position</th><th>Identity</th><th>Conservative similarity</th><th>Second E</th><th>Δ LRR-region start</th><th>LRR location</th><th>Ranking</th></tr></thead><tbody>${[...r.matches].sort((a,b)=>a.start-b.start).map(m=>`<tr><td><button class="motif-link" data-protein="${results.indexOf(r)}" data-hit="${m.start}">${m.sequence}</button></td><td>${m.start}–${m.end}</td><td>${pct(m.identity)}</td><td>${pct(m.conservation)}</td><td>${esc(m.secondE)}</td><td>${m.distanceToFirstLrr??'Unknown'}</td><td>${m.repeatIndices.length?'LRR '+m.repeatIndices.join(', '):m.overlappingRepeatIndices.length?'Overlaps LRR '+m.overlappingRepeatIndices.join(', '):m.inFirstLrr?'LRR 1':m.inLrr===true?(r.domainEvidence.lrrRepeats.length?'Between annotated repeats':'LRR region (repeat unassigned)'):m.inLrr===false?'Outside LRR region':'Unknown'}</td><td>${[m.start===r.bestMatch?.start?'Preferred':'',m.start===r.rankings.mostConserved?.start?'Most conserved':'',m.start===r.rankings.nearest?.start?'Nearest LRR':''].filter(Boolean).join(' · ')}</td></tr>`).join('')||'<tr><td colspan="8">No qualifying hits</td></tr>'}</tbody></table></div>`;}
function interfaceTable(hits){return `<div class="table-wrap"><table><thead><tr><th>Pattern</th><th>Matched sequence</th><th>Position</th><th>Constrained residues</th><th>Context</th></tr></thead><tbody>${hits.map(h=>`<tr><td>${h.pattern}</td><td><code>${h.sequence}</code></td><td>${h.start}–${h.end}</td><td>${h.constrainedPositions.map(p=>`${h.sequence[p-h.start]}${p}`).join(', ')}</td><td>${[h.nearFirstLrr?'Near first LRR':'',h.inLrr&&h.nearNleel?'LRR near preferred NLEEL':''].filter(Boolean).join('; ')||'Elsewhere'}</td></tr>`).join('')||'<tr><td colspan="5">No matches in this region</td></tr>'}</tbody></table></div>`;}
function card(r,index){const m=r.bestMatch,selected=r.matches.find(h=>h.start===selectedHits.get(r.id))??m,c=contextFor(r.sequence,selected),a=annotations.get(r.id)??{},p=phenotypeRecords(r);
 return `<article class="result-card" id="result-${index}"><header class="result-head"><div><p class="eyebrow">${esc(r.id)}</p><h2>${r.length.toLocaleString()} aa · ${esc(r.domainEvidence.nlrClass)}</h2></div>${badge(r.integrated.label)}</header>
 <p class="tiny">Integrated call basis: ${esc(r.integrated.basis)}. Sequence-only call: ${esc(r.classification.label)}.</p>
 <div class="score-strip"><div><span>Preferred candidate · LRR-start priority</span><strong>${esc(loc(m))}</strong><span>${pct(m?.conservation)} conservative similarity · ${pct(m?.identity)} identity</span></div><div><span>Most conserved anywhere</span><strong>${esc(loc(r.rankings.mostConserved))}</strong></div><div><span>Nearest LRR-region start</span><strong>${esc(loc(r.rankings.nearest))}</strong></div></div>
 <section class="interpretation"><h3>Evidence-aware interpretation</h3><p>${esc(r.classification.explanation)}</p><p>${esc(r.integrated.explanation)}</p><p class="tiny">${esc(r.rankings.selectionReason)}. Similarity is not a phenotype probability.</p></section>
 <section><h3>Domain architecture and motif positions</h3><p class="domain-status">${esc(([...domainJobs.values()].filter(j=>j.sequenceId===r.id).map(j=>j.provider+' · '+j.message).join('; ')||null)??(r.domainEvidence.domains.length?'Source annotations matched to this exact sequence':'Domain prediction not run; dashed estimates are not verified domains'))}</p><p><strong>LRR region start: ${r.domainEvidence.lrrRegionStart??'Unresolved'}</strong> · ${esc(r.domainEvidence.lrrRegionSource)}</p>${schematic(r)}</section>
 <details><summary>Set domain coordinates for ${esc(r.id)}</summary><form class="annotation-form" data-index="${index}"><p class="tiny">Use one-based coordinates from a source-supported annotation. WHD is the winged-helix subdomain, not the entire NB-ARC. Unfilled boundaries remain inferred/unavailable.</p><div class="annotation-fields"><label>First LRR start<input name="firstLrrStart" type="number" min="1" max="${r.length}" value="${a.firstLrrStart??''}"></label><label>First LRR end<input name="firstLrrEnd" type="number" min="1" max="${r.length}" value="${a.firstLrrEnd??''}"></label><label>LRR region end<input name="lrrEnd" type="number" min="1" max="${r.length}" value="${a.lrrEnd??''}"></label><label>WHD start<input name="whdStart" type="number" min="1" max="${r.length}" value="${a.whdStart??''}"></label><label>WHD end<input name="whdEnd" type="number" min="1" max="${r.length}" value="${a.whdEnd??''}"></label><label>LRR repeat ranges<textarea name="repeatRanges" placeholder="552-567, 574-589">${esc((a.lrrRepeats??[]).map(x=>`${x.start}-${x.end}`).join(', '))}</textarea></label><label>Annotation source<input name="source" value="${esc(a.source??'')}" placeholder="Publication, domain accession, or structural annotation"></label></div><button class="secondary" type="submit">Apply coordinates</button><span class="annotation-error" role="alert"></span></form></details>
 <section class="interpretation"><h3>LRR region and motif counts</h3><p><strong>Qualifying motif near the annotated LRR-region start: ${r.domainEvidence.boundaryVerified?(r.matches.some(h=>h.proximalVerified)?'Yes':'No'):'Unknown'}</strong> · ±50 aa; individual-repeat containment is reported separately.</p><p><strong>Fully inside individual LRR repeat 1: ${r.rankings.summary.firstLrrStatus}</strong>${r.domainEvidence.firstLrrEnd?` · LRR 1: ${r.domainEvidence.firstLrrStart}–${r.domainEvidence.firstLrrEnd}`:' · Supply first-repeat start and end to establish containment.'}</p><p>${r.rankings.summary.totalHits} total motif hit(s) · ${r.rankings.summary.lrrHits??'Unknown'} in the annotated LRR region · ${r.rankings.summary.firstLrrHits??'Unknown'} inside LRR 1 · ${r.rankings.summary.proximalHits??'Unknown'} near the LRR start (±50 aa) · ${r.rankings.summary.repeatCount??'Unknown'} annotated LRR repeats.</p><p class="tiny">Counts describe qualifying NLEEL-like sequence motifs, not separate protein domains. Near the LRR start is distinct from inside the first repeat. Unknown means missing annotation, not absence.</p></section>
 
 <section id="window-${index}" class="motif-window"><h3>Preferred motif window · ${esc(loc(m))}</h3><p>${esc(r.rankings.selectionReason)}.</p>${m?`<div class="sequence">${numberedSequence(r.sequence.slice(Math.max(0,m.start-51),Math.min(r.length,m.end+50)),[m],Math.max(1,m.start-50))}</div><p class="tiny">Identity ${pct(m.identity)} · conservative similarity ${pct(m.conservation)} · second-E-equivalent ${m.secondEResidue}${m.start+3} · fully in LRR 1: ${m.inFirstLrr==null?'Unknown':m.inFirstLrr?'Yes':'No'}</p>`:'<p>No qualifying motif.</p>'}
 ${selected&&selected.start!==m?.start?`<h3>Selected motif · ${esc(loc(selected))}</h3><div class="sequence">${numberedSequence(r.sequence.slice(c.start-1,c.end),[selected],c.start)}</div><p>Identity ${pct(selected.identity)} · conservative similarity ${pct(selected.conservation)} · second-E-equivalent ${selected.secondEResidue}${selected.start+3} · Δ LRR-region start ${selected.distanceToFirstLrr??'Unknown'}</p>`:''}</section>
 <section><h3>All qualifying NLEEL-like matches</h3>${hitTable(r)}${r.nearStartCandidates.length?`<details><summary>Why nearby windows did not qualify</summary><p class="tiny">The retained rule requires ≥80% conservative similarity and ≥60% identity to the five-residue NLEEL pattern. An EEL substring alone is not sufficient. These near-start windows are shown for inspection and are not counted as qualifying hits.</p><div class="table-wrap"><table><thead><tr><th>Window</th><th>Position</th><th>Identity</th><th>Conservative similarity</th><th>Second-E-equivalent</th></tr></thead><tbody>${r.nearStartCandidates.slice(0,8).map(h=>`<tr><td>${h.sequence}</td><td>${h.start}–${h.end}</td><td>${pct(h.identity)}</td><td>${pct(h.conservation)}</td><td>${h.secondEResidue}${h.start+3}</td></tr>`).join('')}</tbody></table></div></details>`:''}</section>
 <section><h3>Reference-guided positional comparison</h3>${referenceHtml(r.referenceComparisons,r)}</section>
 <section class="interface-box"><h3>WHD–LRR candidate region · sequence evidence</h3><p><code>axxaxxa / LxxL / IxaxxCxaxxaxx</code> · ${r.interfaceEvidence.localHits.length} local segment match(es); ${r.interfaceEvidence.complete?'complete':'incomplete'} exploratory ordered bundle.</p><p class="tiny">a = A/I/L/M/V; x = any standard residue. Slash/order interpretation remains an unverified hypothesis. A hit is not a measured interface.</p>${interfaceTable(r.interfaceEvidence.localHits)}<details><summary>All ${r.interfaceEvidence.hits.length} segment matches and coordinates</summary>${interfaceTable(r.interfaceEvidence.hits)}</details></section>
 <section><h3>Temperature phenotype ledger</h3>${p.length?p.map(x=>`<p><strong>${esc(x.system)} · ${esc(x.direction||'unresolved')}</strong> — ${esc(x.conditions)} ${safeLink(x.reference)}</p><p class="tiny">${esc(x.assay)} · ${esc(x.status)} · ${esc(x.notes)}</p>`).join(''):'<p class="tiny">No exact-sequence curated phenotype. Name or homolog similarity does not transfer a phenotype.</p>'}</section>
 <section><h3>Full numbered sequence · all motif hits highlighted</h3><div class="sequence full">${sequenceHtml(r.sequence,r.matches)}</div></section></article>`;
}
function render(){
 $('#short-reference-check').hidden=!results.length;$('#short-reference-results').innerHTML=results.map(r=>shortReferenceHtml(r)+jointReferenceHtml(r)).join('');
 $('#results').innerHTML=results.map(card).join('');
 $('#batch-body').innerHTML=results.map((r,i)=>`<tr><td><a href="#result-${i}">${esc(r.id)}</a></td><td>${r.length}</td><td>${esc(loc(r.bestMatch))}</td><td>${esc(r.domainEvidence.lrrRegionStart??'Unknown')}</td><td>${r.rankings.summary.totalHits} / ${r.rankings.summary.lrrHits??'Unknown'}</td><td>${r.rankings.summary.firstLrrStatus}</td><td>${r.interfaceEvidence.localHits.length}</td><td>${esc(r.classification.label)}</td><td>${badge(r.integrated.label)}</td></tr>`).join('');
 $('#batch-summary').hidden=!results.length;$('#structure-workspace').hidden=!results.length;
 document.querySelectorAll('[data-hit]').forEach(b=>b.addEventListener('click',()=>{const i=Number(b.dataset.protein);selectedHits.set(results[i].id,Number(b.dataset.hit));structureChoices.set(results[i].id,{choice:b.dataset.hit});render();if(Number($('#structure-sequence').value)===i){syncStructureControls();updateStructureSelection();}document.querySelector(`#window-${i}`).scrollIntoView({block:'center'});}));
 document.querySelectorAll('.annotation-form').forEach(form=>form.addEventListener('submit',e=>{e.preventDefault();try{
   const index=Number(form.dataset.index),old=results[index],a=Object.fromEntries(new FormData(form));
   for(const k of ['firstLrrStart','firstLrrEnd','lrrEnd','whdStart','whdEnd'])if(!a[k])delete a[k];
   if(a.repeatRanges.trim())a.lrrRepeats=a.repeatRanges.split(',').map(x=>{const m=x.trim().match(/^(\d+)\s*[-–]\s*(\d+)$/);if(!m)throw new Error('Enter repeat ranges as start-end, separated by commas.');return {start:Number(m[1]),end:Number(m[2])};});delete a.repeatRanges;
   if(Object.keys(a).some(k=>k!=='source')&&!a.source.trim())throw new Error('Add the annotation source so these coordinates are traceable.');
   a.domains=annotations.get(old.id)?.domains??[];
   for(const k of ['lrrRegionStart','lrrRegionEnd','lrrRegionSource','domainJob'])if(annotations.get(old.id)?.[k]!=null)a[k]=annotations.get(old.id)[k];
   if(annotations.get(old.id)?.nlrClass)a.nlrClass=annotations.get(old.id).nlrClass;
   const next=analyzeRecord(old,config,a);delete next.structure;annotations.set(old.id,a);results[index]=combine(next);
   render();if(Number($('#structure-sequence').value)===index&&payload)updateModel();
 }catch(error){form.querySelector('.annotation-error').textContent=error.message;}}));
}
function run(){
 if(domainBusy)throw new Error('Wait for domain analysis, or stop the queue before starting a new FASTA batch.');
 if(!config)throw new Error('Evidence rules have not loaded.');
 const records=parseFasta($('#fasta-input').value);annotations=new Map();selectedHits.clear();structureChoices.clear();domainJobs.clear();
 results=records.map(record=>{const p=ledger.records.find(p=>p.sequence===record.sequence&&p.annotation);const a=deriveRegionAnnotation(p?.annotation??{});if(p)annotations.set(record.id,a);return combine(analyzeRecord(record,config,a));});
 viewer?.clear();payload=null;chains=[];$('#structure-evidence').textContent='Load a full-length model to assess mapping and contacts.';$('#structure-status').textContent='No structure loaded';$('#residue-context').innerHTML='Load a matching structure to assess motif residues.';
 render();$('#structure-sequence').innerHTML=results.map((r,i)=>`<option value="${i}">${esc(r.id)}</option>`).join('');$('#alphafold-id').value=inferUniProtAccession(results[0]?.description);$('#structure-chain').innerHTML='';syncStructureControls();
 $('#status').textContent=`${results.length} proteins scanned · rules ${config.schemaVersion}`;
 if($('#predict-domains').checked){if($('#domain-provider').value==='cdd'||$('#job-email').value.trim())void predictDomains();else $('#domain-status').textContent='Enter your job email and select Predict domains to run InterProScan for this batch. Motif-only results remain provisional.';}
}
function download(text,name,type='text/plain'){const url=URL.createObjectURL(new Blob([text],{type})),link=document.createElement('a');link.href=url;link.download=name;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function updateModel(){
 const r=results[Number($('#structure-sequence').value)];if(!r||!payload)return;
 if(structureSequence!==r.sequence)throw new Error('Choose a matching model for the selected sequence.');
 const chain=chains[Number($('#structure-chain').value)];if(!chain)throw new Error('Choose a protein chain.');
 const {mapping,measurements}=viewer.display(r,chain,config.structure,$('#confidence-plddt').checked);
 r.structure=measurements;combine(r);render();syncStructureControls();updateStructureSelection();
 $('#structure-status').textContent=`${payload.source} · chain ${chain.id||'(blank)'} · ${mapping.status}`;
 $('#plddt-summary').textContent=$('#confidence-plddt').checked&&measurements.status==='measured'?`Mean pLDDT: WHD ${measurements.whdPlddt?.toFixed(1)??'—'} · LRR ${measurements.lrrPlddt?.toFixed(1)??'—'}`:'Confidence field is not being used as pLDDT, or domain confidence is unavailable.';
 $('#structure-evidence').innerHTML=measurements.status!=='measured'?`<p>${esc(measurements.reason)}</p>${mapping.status==='mapped'?`<p>Sequence coverage: ${pct(mapping.coverage)}. Mapping: ${esc(mapping.method)}.</p>`:''}`:
 `<h3>${measurements.totalPairs} WHD–LRR residue-pair contacts ≤${measurements.cutoff} Å</h3><p>WHD coverage ${pct(measurements.whdCoverage)} · LRR coverage ${pct(measurements.lrrCoverage)}. ${esc(measurements.reason)}</p><div class="table-wrap"><table><thead><tr><th>Local pattern</th><th>Position</th><th>Mapped residues</th><th>Window contact pairs</th><th>Constrained-position pairs</th></tr></thead><tbody>${measurements.localMotifs.map(h=>`<tr><td>${h.pattern}</td><td>${h.start}–${h.end}</td><td>${h.mappedResidues}/${h.sequence.length}</td><td>${h.windowPairs}</td><td>${h.constrainedPairs}</td></tr>`).join('')||'<tr><td colspan="5">No local shorthand segments; domain contact total is shown separately.</td></tr>'}</tbody></table></div><details><summary>Contact positions and minimum heavy-atom distances</summary><div class="table-wrap"><table><thead><tr><th>WHD residue</th><th>LRR residue</th><th>Author coordinates</th><th>Distance Å</th></tr></thead><tbody>${measurements.contacts.map(p=>`<tr><td>${p.whdAA}${p.whdPosition}</td><td>${p.lrrAA}${p.lrrPosition}</td><td>${esc(p.whdAuthor)} ↔ ${esc(p.lrrAuthor)}</td><td>${p.distance.toFixed(2)}</td></tr>`).join('')}</tbody></table></div></details>`;
 $('#download-contacts').disabled=measurements.status!=='measured';
}
async function loadModel(p){viewer??=createStructureViewer($('#viewer'));chains=viewer.load(p.text,p.format);payload=p;structureSequence=results[Number($('#structure-sequence').value)].sequence;$('#confidence-plddt').checked=Boolean(p.confidenceIsPlddt);$('#structure-chain').innerHTML=chains.map((c,i)=>`<option value="${i}">Chain ${esc(c.id||'(blank)')} · ${c.residues.length} aa</option>`).join('');updateModel();}
const structureError=e=>{$('#structure-status').textContent=e.message;};
$('#analysis-form').addEventListener('submit',e=>{e.preventDefault();try{run();}catch(error){$('#status').textContent=error.message;}});
document.querySelectorAll('[data-example]').forEach(b=>b.addEventListener('click',()=>{$('#fasta-input').value=EXAMPLES[b.dataset.example].fasta;}));
$('#fasta-file').addEventListener('change',async e=>{try{const f=e.target.files[0];if(f){if(f.size>10e6)throw new Error('Limit FASTA uploads to 10 MB.');$('#fasta-input').value=await f.text();$('#status').textContent=`${f.name} loaded`;}}catch(error){$('#status').textContent=error.message;}});
$('#download-csv').addEventListener('click',()=>download(resultsTable(results,','),'nlr-temper-results.csv','text/csv'));
$('#download-tsv').addEventListener('click',()=>download(resultsTable(results,'\t'),'nlr-temper-results.tsv','text/tab-separated-values'));
$('#structure-sequence').addEventListener('change',()=>{viewer?.clear();payload=null;chains=[];$('#structure-chain').innerHTML='';$('#alphafold-id').value=inferUniProtAccession(results[Number($('#structure-sequence').value)].description);$('#structure-status').textContent='Load a model for this sequence';$('#structure-evidence').textContent='No model loaded for the selected sequence.';$('#residue-context').innerHTML='Load a matching structure to assess motif residues.';$('#plddt-summary').textContent='Confidence unavailable';$('#download-contacts').disabled=true;syncStructureControls();});
$('#structure-chain').addEventListener('change',()=>{try{updateModel();}catch(e){structureError(e);}});
$('#confidence-plddt').addEventListener('change',()=>{try{updateModel();}catch(e){structureError(e);}});
$('#load-alphafold').addEventListener('click',async()=>{try{$('#load-alphafold').disabled=true;$('#structure-status').textContent='Retrieving and checking AlphaFold sequence…';await loadModel(await fetchAlphaFoldPdb($('#alphafold-id').value,results[Number($('#structure-sequence').value)].sequence));}catch(e){structureError(e);}finally{$('#load-alphafold').disabled=false;}});
$('#structure-file').addEventListener('change',async e=>{try{const f=e.target.files[0];if(!f)return;if(f.size>40e6)throw new Error('Limit structure files to 40 MB.');await loadModel({text:await f.text(),format:/\.(pdb|ent)$/i.test(f.name)?'pdb':'cif',source:f.name,confidenceIsPlddt:false});}catch(e){structureError(e);}});
$('#focus-motif').addEventListener('click',()=>viewer?.focus(results[Number($('#structure-sequence').value)]));$('#reset-view').addEventListener('click',()=>viewer?.reset());
$('#download-contacts').addEventListener('click',()=>{const r=results[Number($('#structure-sequence').value)],s=r.structure;if(s?.status!=='measured')return;download(delimited([['sequence_id','chain','whd_position','whd_aa','lrr_position','lrr_aa','minimum_heavy_atom_distance_A','whd_author','lrr_author','cutoff_A'],...s.contacts.map(p=>[r.id,s.chain,p.whdPosition,p.whdAA,p.lrrPosition,p.lrrAA,p.distance,p.whdAuthor,p.lrrAuthor,s.cutoff])]),'nlr-temper-contacts.tsv');});
window.addEventListener('resize',()=>viewer?.resize());
async function init(){try{
 const [rules,phenotypes,validation,referenceData,domainData]=await Promise.all([fetch('./config/motifs.v0.5.0.json'),fetch('./data/phenotypes.json'),fetch('./data/validation-report.json'),fetch('./data/position-references.json'),fetch('./data/reference-domains.json')]);
 if(!rules.ok)throw new Error('Evidence configuration could not load.');config=await rules.json();
 if(referenceData.ok)references=await referenceData.json();
 if(phenotypes.ok)ledger=await phenotypes.json();if(validation.ok)report=await validation.json();
 if(domainData.ok){const domains=await domainData.json();for(const r of ledger.records)if(r.annotation&&domains[r.system])r.annotation.domains=domains[r.system];}
 $('#analyze-button').disabled=false;$('#status').textContent=`Frozen rules ${config.schemaVersion} · local analysis`;
 if(report)$('#validation-summary').innerHTML=`<p>${esc(report.summary)}</p><p>${esc(report.validationStatus)}</p><div class="download-actions"><a href="./data/validation-report.json" download>Validation report (JSON)</a><a href="./data/phenotypes.json" download>Phenotype ledger</a><a href="./data/benchmark-results.tsv" download>Benchmark results</a><a href="./data/splits.tsv" download>Homology groups / holdouts</a></div>`;
 }catch(e){$('#status').textContent=e.message;}}
init();

$('#domain-file').addEventListener('change',async e=>{try{if(!results.length)throw new Error('Analyze FASTA sequences before importing their annotations.');const file=e.target.files[0];if(!file)return;if(file.size>20e6)throw new Error('Limit annotations to 20 MB.');const imported=parseDomains(await file.text(),$('#domain-format').value,results),nextAnnotations=new Map(annotations);const next=results.map(r=>{if(!imported.has(r.id))return r;const a={...(annotations.get(r.id)??{}),...imported.get(r.id)};nextAnnotations.set(r.id,a);const n=analyzeRecord(r,config,a);delete n.structure;return combine(n);});annotations=nextAnnotations;results=next;render();if(payload)updateModel();$('#domain-status').textContent=`Imported source annotations for ${imported.size} protein(s). Generic LRR signature spans do not establish individual repeats; use exact-sequence JSON or the coordinate form for repeat annotation.`;}catch(e){$('#domain-status').textContent=e.message;}finally{e.target.value='';}});

function showDomainJobs(){
 $('#domain-progress').innerHTML=[...domainJobs.values()].map(j=>`<p><strong>${esc(j.sequenceId)} · ${esc(j.provider)}</strong> · ${esc(j.message)}${j.id?` <code>${esc(j.id)}</code>`:''}</p>`).join('');
}
async function predictDomains(){
 if(domainBusy)return;if(!results.length){$('#domain-status').textContent='Analyze FASTA first.';return;}
 const choice=$('#domain-provider').value||'interpro',providers=choice==='both'?['cdd','interpro']:[choice],email=$('#job-email').value.trim();
 if(providers.includes('interpro')&&!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)){$('#domain-status').textContent='Enter a job email for InterProScan, or select NCBI CDD to run without email.';return;}
 domainBusy=true;stopDomainJobs=false;$('#predict-domain-button').disabled=true;$('#stop-domain-button').disabled=false;$('#analyze-button').disabled=true;$('#domain-provider').disabled=true;
 $('#domain-status').textContent='Running selected domain services. Short reference alignments are already available above. Keep this tab open.';
 try{for(const provider of providers){for(let index=0;index<results.length;index++){
  if(stopDomainJobs)break;const r=results[index],key=provider+':'+r.id;let job=domainJobs.get(key);
  if(job?.done)continue;
  try{
   if(!job){job={sequenceId:r.id,provider:provider==='cdd'?'NCBI CDD':'InterProScan',message:'Submitting…'};domainJobs.set(key,job);}
   if(!job.id){showDomainJobs();if(provider==='cdd'){job.response=await cddRequest({action:'submit',sequence:r.sequence});job.id=job.response.id;if(!job.id)throw new Error('CDD did not return a job ID.');}else job.id=await submitDomainJob(r,email);}
   let status='PENDING',polls=0;while(!stopDomainJobs){
    if(provider==='cdd'){job.response=job.response?.status==='FINISHED'?job.response:await cddRequest({action:'poll',id:job.id});status=job.response.status;}else status=await domainJobStatus(job.id);
    job.message=status==='FINISHED'?'Collecting domain annotations':status;showDomainJobs();if(status==='FINISHED')break;
    if(!['RUNNING','PENDING','QUEUED'].includes(status))throw new Error(`Service status: ${status}.`);
    if(++polls>=360)throw new Error('Still queued after one hour. Select Predict domains later to resume.');await new Promise(resolve=>setTimeout(resolve,10000));
   }
   if(stopDomainJobs){job.message='Monitoring paused; remote job may continue. Select Predict domains to resume.';break;}
   const prediction=provider==='cdd'?collectCddResult(job.response.text,r,job.id):await collectDomainJob(job.id,r);
   prediction.domains=prediction.domains.map(d=>({...d,engine:provider}));const previous=annotations.get(r.id)??{};
   const a=deriveRegionAnnotation({...previous,...prediction,domains:[...(previous.domains??[]).filter(d=>d.engine!==provider),...prediction.domains]});
   annotations.set(r.id,a);const next=analyzeRecord(r,config,a);delete next.structure;delete next.residueContext;results[index]=combine(next);
   job.done=true;job.message=prediction.domains.length?`Complete · ${prediction.domains.length} signature matches; annotations applied`:'Complete · no signature matches reported';render();showDomainJobs();
   if(payload&&Number($('#structure-sequence').value)===index)updateModel();
  }catch(error){job??={sequenceId:r.id,provider};job.message=error.message;domainJobs.set(key,job);showDomainJobs();break;}
 }if(stopDomainJobs)break;}}finally{domainBusy=false;$('#predict-domain-button').disabled=false;$('#stop-domain-button').disabled=true;$('#analyze-button').disabled=false;$('#domain-provider').disabled=false;$('#domain-status').textContent=stopDomainJobs?'Queue paused. Existing results retained.':'Domain queue stopped or complete. Check each service status below; unfinished jobs can be resumed.';showDomainJobs();}
}
$('#predict-domain-button').addEventListener('click',predictDomains);
$('#stop-domain-button').addEventListener('click',()=>{stopDomainJobs=true;$('#domain-status').textContent='Pausing after the current request. Submitted remote jobs are not cancelled.';});
function renderResidueContext(r){const context=r.residueContext;
 $('#residue-context').innerHTML=context?.status==='measured'?`<h3>Selected residues: secondary structure and exposure</h3><p class="tiny">${esc(context.method)}</p><div class="table-wrap"><table><thead><tr><th>Residue</th><th>Motif</th><th>Secondary structure</th><th>Outside coil/turn?</th><th>Approx. SASA Å²</th><th>Heavy atoms</th><th>pLDDT</th><th>WHD–LRR pairs</th></tr></thead><tbody>${context.rows.map(x=>`<tr><td><button class="motif-link" data-residue="${x.position}">${x.aa}${x.position}${r.matches.some(m=>m.start+3===x.position)?' · second-E equivalent':''}</button></td><td>${esc(x.motifs.join('; '))}</td><td>${esc(x.secondary)}</td><td>${x.outsideCoil??'Unknown'}</td><td>${x.sasa?.toFixed(1)??'Unknown'}</td><td>${x.atomsComplete?'Complete count':'Incomplete / unknown'}</td><td>${x.confidence?.toFixed(1)??'Unavailable'}${x.confidence!=null&&x.confidence<70?' · low confidence':''}</td><td>${x.contacts??'Unresolved'}</td></tr>`).join('')}</tbody></table></div>`:`<p>${esc(context?.reason??'Load a matching model.')}</p>`;
 document.querySelectorAll('[data-residue]').forEach(b=>b.addEventListener('click',()=>viewer?.focusResidue(Number(b.dataset.residue))));
}
$('#download-residues').addEventListener('click',()=>{const r=results[Number($('#structure-sequence').value)];if(!r?.residueContext?.rows.length)return;download(delimited([['sequence','position','aa','motifs','secondary_structure','outside_coil','approx_SASA_A2','complete_heavy_atom_count','plddt','whd_lrr_pairs'],...r.residueContext.rows.map(x=>[r.id,x.position,x.aa,x.motifs.join(';'),x.secondary,x.outsideCoil,x.sasa,x.atomsComplete,x.confidence,x.contacts])]),'nlr-temper-motif-structure.tsv');});

function syncStructureControls(){
 const r=results[Number($('#structure-sequence').value)];if(!r)return;
 let state=structureChoices.get(r.id)??{choice:'preferred'};
 if(!['preferred','all','custom'].includes(state.choice)&&!r.matches.some(m=>String(m.start)===state.choice))state={choice:'preferred'};
 structureChoices.set(r.id,state);
 $('#structure-motif').innerHTML=`<option value="preferred">Automatic preferred · ${esc(loc(r.bestMatch))}</option><option value="all">All ${r.matches.length} detected motifs</option>${r.matches.map(m=>`<option value="${m.start}">${esc(loc(m))}${m.start===r.bestMatch?.start?' · preferred':''}</option>`).join('')}<option value="custom">Custom residue range</option>`;
 $('#structure-motif').value=state.choice;
 $('#custom-region').hidden=state.choice!=='custom';
 $('#region-start').max=r.length;$('#region-end').max=r.length;
 $('#region-start').value=state.start??r.bestMatch?.start??1;$('#region-end').value=state.end??r.bestMatch?.end??Math.min(5,r.length);
 if(!payload){$('#selection-summary').textContent='Choose a motif or range, then load a matching structure.';$('#selection-contacts').innerHTML='';}
}
function updateStructureSelection(){
 const r=results[Number($('#structure-sequence').value)];if(!r)return;
 const state=structureChoices.get(r.id)??{choice:'preferred'};
 try{
  const selection=resolveStructureSelection(r,state.choice,state.start??$('#region-start').value,state.end??$('#region-end').value);
  if(!payload||!viewer){$('#selection-summary').textContent=`${selection.label}. Load a matching model to highlight these residues.`;return;}
  const {mapping,nearby}=viewer.highlight(selection,$('#residue-labels').value||'off');
  $('#selection-summary').textContent=`${selection.label} · ${nearby.status==='measured'?`${nearby.mapped}/${nearby.total} residues mapped`:nearby.reason}. Inspection selection; automatic ranking is unchanged.`;
  r.residueContext=residueContext(mapping,{...r,matches:selection.ranges.map(m=>m.custom?{...m,sequence:'Custom region'}:m)},{confidenceIsPlddt:$('#confidence-plddt').checked,measurements:r.structure});
  renderResidueContext(r);renderSelectionContacts(nearby);
 }catch(e){$('#selection-summary').textContent=e.message;}
}
function renderSelectionContacts(nearby){
 if(nearby.status!=='measured'){$('#selection-contacts').innerHTML='';return;}
 const residue=r=>r.position==null?`${r.aa} · ${esc(authorCoordinate(r))}`:`<button class="motif-link" data-neighbor-position="${r.position}">${r.aa}${r.position}</button>`;
 $('#selection-contacts').innerHTML=`<h3>Residues near the selected region</h3><p>${nearby.contacts.length} residue-pair contact(s) with ${nearby.neighbors.length} other protein residue(s) within ${nearby.cutoff} Å.</p><p class="tiny">Minimum heavy-atom distances in the loaded model, across all loaded protein chains. Residues inside the selected region are excluded as neighbors; adjacent sequence residues are included and identified by separation. Other chains retain their PDB numbering. This geometric view is independent of the annotated WHD–LRR contact analysis below.</p>${nearby.contacts.length?`<details><summary>Inspect contact positions and distances</summary><div class="table-wrap"><table><thead><tr><th>Selected residue</th><th>Neighbor</th><th>PDB coordinates</th><th>Sequence separation</th><th>Distance Å</th></tr></thead><tbody>${nearby.contacts.map(p=>`<tr><td>${residue(p.selected)}</td><td>${residue(p.neighbor)}</td><td>${esc(authorCoordinate(p.selected))} ↔ ${esc(authorCoordinate(p.neighbor))}</td><td>${p.separation??'Other chain'}</td><td>${p.distance.toFixed(2)}</td></tr>`).join('')}</tbody></table></div></details>`:'<p>No neighboring residue contacts detected in the mapped coordinates at this cutoff.</p>'}`;
 document.querySelectorAll('[data-neighbor-position]').forEach(b=>b.addEventListener('click',()=>viewer?.focusResidue(Number(b.dataset.neighborPosition))));
}
$('#structure-motif').addEventListener('change',()=>{
 const r=results[Number($('#structure-sequence').value)];if(!r)return;
 const choice=$('#structure-motif').value;structureChoices.set(r.id,{choice,start:Number($('#region-start').value),end:Number($('#region-end').value)});
 if(!['preferred','all','custom'].includes(choice))selectedHits.set(r.id,Number(choice));else if(choice==='preferred')selectedHits.delete(r.id);
 syncStructureControls();updateStructureSelection();render();
});
for(const id of ['#region-start','#region-end'])$(id).addEventListener('input',()=>{const r=results[Number($('#structure-sequence').value)];if(!r)return;structureChoices.set(r.id,{choice:'custom',start:Number($('#region-start').value),end:Number($('#region-end').value)});updateStructureSelection();});
$('#residue-labels').addEventListener('change',()=>viewer?.setLabels($('#residue-labels').value));
