const STANDARD = "ACDEFGHIKLMNPQRSTVWY";
const AMINO_ACIDS = new Set(`${STANDARD}XBZJUO`);

export function parseFasta(text, { maxRecords = 500 } = {}) {
  const normalized = String(text ?? "").replace(/\r/g, "").trim();
  if (!normalized) throw new Error("Paste or upload a protein sequence.");
  const records = []; let current = null;
  for (const raw of normalized.split("\n")) {
    const line = raw.trim(); if (!line || line.startsWith(";")) continue;
    if (line.startsWith(">")) {
      if (current) records.push(current);
      const description = line.slice(1).trim();
      if (!description) throw new Error("A FASTA header needs an identifier.");
      current = { id: description.split(/\s+/)[0], description, sequence: "" };
    } else {
      current ??= { id: "sequence_1", description: "sequence_1", sequence: "" };
      current.sequence += line.replace(/\s+/g, "").toUpperCase();
    }
  }
  if (current) records.push(current);
  if (records.length > maxRecords) throw new Error(`A browser run accepts up to ${maxRecords} records.`);
  if (records.reduce((n,r) => n+r.sequence.length,0) > 3000000 && maxRecords !== Infinity) throw new Error("Limit each browser run to 3 million residues.");
  const seen = new Set();
  return records.map(record => {
    record.sequence = record.sequence.replace(/\*$/, "");
    if (!record.sequence) throw new Error(`${record.id} has no sequence.`);
    const invalid = [...new Set([...record.sequence].filter(aa => !AMINO_ACIDS.has(aa)))];
    if (invalid.length) throw new Error(`${record.id} contains unsupported residue(s): ${invalid.join(", ")}.`);
    if (record.sequence.length > 15000 && maxRecords !== Infinity) throw new Error(`${record.id} exceeds the 15,000-aa interactive limit.`);
    let id = record.id, suffix = 2; while (seen.has(id)) id = `${record.id}_${suffix++}`; seen.add(id);
    return { ...record, id, ambiguousResidues: [...record.sequence].filter(aa => !STANDARD.includes(aa)).length };
  });
}

export function scoreWindow(window, definition) {
  const identity = [...window].filter((aa,i) => aa === definition.reference[i]).length / definition.reference.length;
  const conservation = [...window].filter((aa,i) => definition.conservativeGroups[i]?.includes(aa)).length / definition.reference.length;
  const secondEResidue = window[definition.secondEIndex];
  return { sequence:window, identity, conservation, secondEResidue, secondE:secondEResidue === "E" ? "retained (E)" : secondEResidue === "D" ? "acidic substitution (D)" : `not retained (${secondEResidue})` };
}
export function scanNleel(sequence, definition) {
  const matches = [], width = definition.reference.length;
  for (let i=0;i<=sequence.length-width;i++) {
    const score = scoreWindow(sequence.slice(i,i+width),definition);
    if (score.identity >= definition.minimumIdentity && score.conservation >= definition.minimumConservation)
      matches.push({ ...score, start:i+1, end:i+width, exact:score.sequence === definition.reference });
  }
  return matches.sort((a,b) => b.conservation-a.conservation || b.identity-a.identity || a.start-b.start);
}

export function inferDomainEvidence(sequence) {
  const pLoop = /G....GK[ST]/.exec(sequence);
  const after = (re, offset) => { const m=re.exec(sequence.slice(offset)); return m ? { sequence:m[0],start:offset+m.index+1,end:offset+m.index+m[0].length } : null; };
  const p = pLoop ? {sequence:pLoop[0],start:pLoop.index+1,end:pLoop.index+pLoop[0].length}:null;
  const glpl = p ? after(/G[LIVM]P[LIVM]/,p.end) : null;
  const mhd = glpl ? after(/M[HQR][DN]/,glpl.end) : null;
  // Anchors establish ordering only. They are never presented as verified domains.
  const firstLrrStart = mhd && mhd.end+31 <= sequence.length ? mhd.end+31 : null;
  const nbArc = p && mhd ? {start:Math.max(1,p.start-35),end:mhd.end+30}:null;
  return {level:p && glpl ? "supported" : p ? "partial" : "insufficient",anchors:{pLoop:p,glpl,mhd},
    firstLrrStart,boundarySource:firstLrrStart ? "inferred: ordered P-loop/GLPL/MHD + 30 aa; not first-repeat annotation" : "unavailable",
    boundaryVerified:false,lrrEndAnnotated:false,whd:null,nlrClass:"Not inferred",
    heuristicDomains:{nTerminal:nbArc?{start:1,end:nbArc.start-1}:null,nbArc,lrr:firstLrrStart?{start:firstLrrStart,end:sequence.length}:null}};
}

export function applyDomainAnnotation(sequence, inferred, annotation = {}) {
  annotation = {...annotation};
  const result = structuredClone(inferred);
  const coordinate = (value,name) => { const n=Number(value); if (!Number.isInteger(n)||n<1||n>sequence.length) throw new Error(`${name} must be a one-based coordinate within this protein.`);return n; };
  result.lrrRepeats=[]; result.firstLrrEnd=null;
  result.domains=(annotation.domains??[]).map(d=>{const start=coordinate(d.start,'Domain start'),end=coordinate(d.end,'Domain end');if(end<start)throw new Error('Domain end precedes start.');return {...d,start,end};});
  if(annotation.lrrRepeats?.length){
    result.lrrRepeats=annotation.lrrRepeats.map(r=>({start:coordinate(r.start,'LRR repeat start'),end:coordinate(r.end,'LRR repeat end')})).sort((a,b)=>a.start-b.start||a.end-b.end);
    if(result.lrrRepeats.some(r=>r.end<r.start))throw new Error('LRR repeat end precedes its start.');
    const first=result.lrrRepeats[0],last=result.lrrRepeats.at(-1);
    for(const [key,value] of [['firstLrrStart',first.start],['firstLrrEnd',first.end],['lrrEnd',Math.max(...result.lrrRepeats.map(r=>r.end))]]){
      if(annotation[key]!=null&&annotation[key]!==''&&Number(annotation[key])!==value)throw new Error(`${key} conflicts with supplied repeat ranges.`);
      annotation[key]=value;
    }
  }
  if (annotation.firstLrrStart != null && annotation.firstLrrStart !== "") {
    result.firstLrrStart=coordinate(annotation.firstLrrStart,"First LRR start");
    result.boundarySource=annotation.source || "user-supplied annotation";
    result.boundaryVerified=true;
    const end=annotation.lrrEnd ? coordinate(annotation.lrrEnd,"LRR end"):sequence.length;
    if(end<result.firstLrrStart)throw new Error("LRR end precedes its start.");
    result.heuristicDomains.lrr={start:result.firstLrrStart,end};
    result.lrrEndAnnotated=Boolean(annotation.lrrEnd);
  }
  if(annotation.firstLrrEnd!=null&&annotation.firstLrrEnd!==''){
    if(!result.boundaryVerified)throw new Error('First LRR end requires an annotated first LRR start.');
    result.firstLrrEnd=coordinate(annotation.firstLrrEnd,'First LRR end');
    if(result.firstLrrEnd<result.firstLrrStart||result.firstLrrEnd>result.heuristicDomains.lrr.end)throw new Error('First LRR end must lie within the annotated LRR region.');
  }
  if (annotation.whdStart || annotation.whdEnd) {
    const start=coordinate(annotation.whdStart,"WHD start"),end=coordinate(annotation.whdEnd,"WHD end");
    if(end<start)throw new Error("WHD end precedes its start.");
    if(result.firstLrrStart && end>=result.firstLrrStart)throw new Error("WHD and LRR ranges must not overlap.");
    result.whd={start,end,source:annotation.source || "user-supplied annotation"};
  }
  result.repeatBoundaryVerified=result.boundaryVerified;
  result.repeatBoundarySource=result.boundarySource;
  result.lrrRegionStart=result.boundaryVerified?result.firstLrrStart:null;
  result.lrrRegionEnd=result.lrrEndAnnotated?result.heuristicDomains.lrr.end:null;
  result.lrrRegionSource=result.boundarySource;
  if(annotation.lrrRegionStart!=null){
    result.lrrRegionStart=coordinate(annotation.lrrRegionStart,'LRR region start');
    result.lrrRegionEnd=coordinate(annotation.lrrRegionEnd,'LRR region end');
    if(result.lrrRegionEnd<result.lrrRegionStart)throw new Error('LRR region end precedes start.');
    result.lrrRegionSource=annotation.lrrRegionSource||annotation.source||'domain prediction';
    result.boundaryVerified=true;result.lrrEndAnnotated=true;result.boundarySource=result.lrrRegionSource;
    result.heuristicDomains.lrr={start:result.lrrRegionStart,end:result.lrrRegionEnd};
  }
  result.domainJob=annotation.domainJob??null;
  if(annotation.nlrClass)result.nlrClass=annotation.nlrClass;
  return result;
}

export function rankNleel(matches, firstLrrStart, flank=50, domains={}) {
  const firstKnown=Boolean((domains.repeatBoundaryVerified??domains.boundaryVerified)&&domains.firstLrrEnd!=null);
  const regionStart=domains.lrrRegionStart??firstLrrStart;
  const regionKnown=Boolean(domains.boundaryVerified&&domains.lrrEndAnnotated);
  const repeats=domains.lrrRepeats??[];
  const all=matches.map(m=>({...m,distanceToFirstLrr:regionStart==null?null:m.start-regionStart,
    inFirstLrr:firstKnown?m.start>=firstLrrStart&&m.end<=domains.firstLrrEnd:null,
    overlapsFirstLrr:firstKnown?m.end>=firstLrrStart&&m.start<=domains.firstLrrEnd:null,
    inLrr:regionKnown?m.start>=regionStart&&m.end<=domains.heuristicDomains.lrr.end:null,
    repeatIndices:repeats.flatMap((r,i)=>m.start>=r.start&&m.end<=r.end?[i+1]:[]),
    overlappingRepeatIndices:repeats.flatMap((r,i)=>m.end>=r.start&&m.start<=r.end?[i+1]:[]),
    proximal:regionStart==null?null:m.end>=regionStart-flank && m.start<=regionStart+flank,
    proximalVerified:Boolean(domains.boundaryVerified&&regionStart!=null&&m.end>=regionStart-flank&&m.start<=regionStart+flank)}));
  const byPosition=[...all].sort((a,b)=>a.start-b.start);
  const nearest=regionStart==null?null:[...all].sort((a,b)=>Math.abs(a.distanceToFirstLrr)-Math.abs(b.distanceToFirstLrr)||b.conservation-a.conservation||a.start-b.start)[0]??null;
  const firstProximal=byPosition.find(m=>m.proximal)??null;
  const firstLrrHits=byPosition.filter(m=>m.inFirstLrr),preferredFirstLrr=firstLrrHits[0]??null;
  const firstInRegion=byPosition.find(m=>m.inLrr)??null;
  const domainStartPolicy=domains.priority==='domain-start';
  const preferred=domainStartPolicy?(byPosition.find(m=>m.proximalVerified)??firstInRegion??all[0]??null):(preferredFirstLrr??firstInRegion??firstProximal??all[0]??null);
  return {all,mostConserved:all[0]??null,nearest,firstProximal,preferredFirstLrr,firstInRegion,preferred,
    summary:{totalHits:all.length,lrrHits:regionKnown?all.filter(m=>m.inLrr).length:null,repeatCount:repeats.length||null,firstLrrHits:firstKnown?firstLrrHits.length:null,firstLrrStatus:firstKnown?(firstLrrHits.length?'Yes':'No'):'Unknown',proximalHits:regionStart==null?null:all.filter(m=>m.proximal).length},
    selectionReason:domainStartPolicy?(preferred?.proximalVerified?`Earliest qualifying ≥80% conservative-similarity motif within ±${flank} aa of the annotated LRR region start; individual-repeat containment is not required`:firstInRegion?'Earliest qualifying motif in the annotated LRR region; no qualifying near-start motif':'Most conserved qualifying motif; annotated LRR-start support is unresolved or absent'):preferredFirstLrr?'First priority: earliest qualifying ≥80% conservative-similarity hit fully inside annotated LRR 1':firstInRegion?'Earliest qualifying hit fully inside the annotated LRR region; no qualifying LRR 1 hit':firstProximal?`Earliest qualifying hit in the first-LRR ±${flank}-aa window; not confirmed inside LRR 1`:firstLrrStart==null?"Most conserved hit; first-LRR boundary unavailable":"Most conserved hit; no qualifying first-repeat or proximal hit"};
}

export function scanInterfaceSignature(sequence, definition, lrrStart=1, nleel=null, flank=50) {
  const start=lrrStart??1;
  const hitSets=definition.segments.map(pattern=>{
    const hits=[];
    for(let i=0;i<=sequence.length-pattern.length;i++){
      const window=sequence.slice(i,i+pattern.length);
      if(![...window].every(aa=>STANDARD.includes(aa)))continue;
      if(![...pattern].every((p,j)=>p==="x"||(p==="a"?definition.aliphaticResidues.includes(window[j]):p===window[j])))continue;
      hits.push({pattern,sequence:window,start:i+1,end:i+pattern.length,
        constrainedPositions:[...pattern].flatMap((p,j)=>p==="x"?[]:[i+j+1]),
        inLrr:lrrStart!=null && i+1>=lrrStart,
        nearFirstLrr:lrrStart!=null && i+pattern.length>=lrrStart-flank && i+1<=lrrStart+flank,
        nearNleel:Boolean(nleel && i+pattern.length>=nleel.start-flank && i+1<=nleel.end+flank)});
    }return hits;
  });
  let bundle=null;
  // Only the nearest subsequent hits are needed for a minimal-span ordered bundle.
  for(const first of hitSets[0].filter(h=>h.start>=start)) {
    const second=hitSets[1].find(h=>h.start>first.end);if(!second)continue;
    const third=hitSets[2].find(h=>h.start>second.end);if(!third)continue;
    const span=third.end-first.start+1;
    if(span<=definition.maximumBundleSpan && (!bundle||span<bundle.span))bundle={hits:[first,second,third],span};
  }
  const hits=hitSets.flat().sort((a,b)=>a.start-b.start);
  return {complete:Boolean(bundle),score:bundle?1:hitSets.filter(h=>h.some(x=>x.start>=start)).length/definition.segments.length,
    hits,localHits:hits.filter(h=>h.nearFirstLrr || (h.inLrr&&h.nearNleel)),bundle,span:bundle?.span??null,searchStart:start};
}

export function nleelEvidenceScore(match) {
  if(!match)return 0;
  return (0.35*match.identity+0.65*match.conservation)*(["E","D"].includes(match.secondEResidue)?1:0.6);
}
export function classifyEvidence({bestMatch,interfaceEvidence,domainEvidence,config}) {
  const nleelScore=nleelEvidenceScore(bestMatch),interfaceScore=interfaceEvidence.score;
  const base={confidence:"Low — exploratory sequence evidence",nleelScore,interfaceScore};
  if(domainEvidence.level==="insufficient" && !domainEvidence.boundaryVerified)return {...base,label:"Uncertain",explanation:"NLR context is not established; motif matches are not classified."};
  if(!(bestMatch?.inFirstLrr||bestMatch?.proximal) && nleelScore>=config.classification.tsLikeMinimum)return {...base,label:"Uncertain",explanation:"NLEEL-like similarity is strong, but a proximal first-LRR location is unestablished. Whole-protein similarity alone is insufficient."};
  const positionalSupport=config.position?.priority==='domain-start'?bestMatch?.proximalVerified:(bestMatch?.inFirstLrr||bestMatch?.proximal);
  const ts=positionalSupport && nleelScore>=config.classification.tsLikeMinimum;
  if(ts && interfaceEvidence.complete)return {...base,label:"Uncertain",explanation:"Proximal NLEEL-like support competes with the exploratory ordered interface bundle. Geometry and phenotype evidence require separate assessment."};
  if(ts)return {...base,label:"TS-like",explanation:`${bestMatch.sequence} at ${bestMatch.start}–${bestMatch.end} is ${config.position?.priority==='domain-start'?'near the annotated LRR domain-region start':bestMatch.inFirstLrr?'inside annotated LRR 1':`near the ${domainEvidence.boundaryVerified?'supplied':'inferred'} LRR start; first-repeat containment is not established`}; similarity support ${Math.round(nleelScore*100)}%. This is a sequence hypothesis.`};
  return {...base,label:"Uncertain",explanation:interfaceEvidence.complete?"The exploratory interface bundle is present. Sequence presence alone cannot establish WHD–LRR stabilization or temperature tolerance.":`A qualifying sequence match is ${bestMatch?`${bestMatch.sequence} at ${bestMatch.start}–${bestMatch.end}; conservative similarity ${Math.round(bestMatch.conservation*100)}%, identity ${Math.round(bestMatch.identity*100)}%, offset from LRR-region start ${bestMatch.distanceToFirstLrr??'unknown'} aa`:"absent"}. The retained TS-like evidence threshold and annotated near-start support are not both satisfied. Motif absence does not imply temperature tolerance.`};
}
export function integratedCall(result, phenotypes=[], structure=null) {
  const trusted=phenotypes.filter(p=>p.sequenceVerified && p.status==="curated" && p.direction);
  const directions=[...new Set(trusted.map(p=>p.direction))];
  const structuralText=structure?.status==="measured"?`Mapped WHD–LRR geometry: ${structure.totalPairs} residue pairs; this has no calibrated temperature-response direction.`:`WHD–LRR structure evidence: ${structure?.reason||"not assessed"}.`;
  let label=result.classification.label, basis="computational hypothesis";
  if(directions.length===1){label=directions[0]==="heat-sensitive"?"TS-like":directions[0]==="heat-tolerant"||directions[0]==="heat-enhanced"?"TT-like":"Uncertain";basis="published phenotype in stated conditions";}
  if(directions.length>1){label="Uncertain";basis="context-dependent or conflicting phenotype records";}
  return {label,basis,explanation:`Sequence call: ${result.classification.label}. ${structuralText} ${trusted.length?`${trusted.length} exact-sequence phenotype record(s); retain assay, host and temperature context.`:"No curated exact-sequence temperature phenotype. No phenotype is transferred by name or homology."}`,phenotypes:trusted,structure};
}
export function analyzeRecord(record,config,annotation={}) {
  const domainEvidence=applyDomainAnnotation(record.sequence,inferDomainEvidence(record.sequence),annotation);
  domainEvidence.priority=config.position?.priority??'first-repeat';
  const rankings=rankNleel(scanNleel(record.sequence,config.nleel),domainEvidence.firstLrrStart,config.position?.contextFlank??50,domainEvidence);
  const bestMatch=rankings.preferred;
  const interfaceEvidence=scanInterfaceSignature(record.sequence,config.interface,domainEvidence.lrrRegionStart??domainEvidence.firstLrrStart,bestMatch,config.position?.contextFlank??50);
  const nearStartCandidates=[];
  const start=domainEvidence.lrrRegionStart;
  if(start!=null)for(let i=Math.max(0,start-51);i<=Math.min(record.sequence.length-5,start+49);i++){
    const hit=scoreWindow(record.sequence.slice(i,i+5),config.nleel);
    if(hit.conservation>=.6&&!rankings.all.some(m=>m.start===i+1))nearStartCandidates.push({...hit,start:i+1,end:i+5});
  }
  nearStartCandidates.sort((a,b)=>b.conservation-a.conservation||b.identity-a.identity||a.start-b.start);
  return {...record,nearStartCandidates,length:record.sequence.length,ruleset:config.schemaVersion,matches:rankings.all,rankings,bestMatch,domainEvidence,interfaceEvidence,
    classification:classifyEvidence({bestMatch,interfaceEvidence,domainEvidence,config})};
}
export function contextFor(sequence,match,flank=50) {
  if(!match)return null;
  const start=Math.max(1,match.start-flank),end=Math.min(sequence.length,match.end+flank);
  return {before:sequence.slice(start-1,match.start-1),motif:sequence.slice(match.start-1,match.end),after:sequence.slice(match.end,end),start,end};
}
