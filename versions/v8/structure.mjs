const AA={ALA:"A",ARG:"R",ASN:"N",ASP:"D",CYS:"C",GLN:"Q",GLU:"E",GLY:"G",HIS:"H",ILE:"I",LEU:"L",LYS:"K",MET:"M",PHE:"F",PRO:"P",SER:"S",THR:"T",TRP:"W",TYR:"Y",VAL:"V",MSE:"M"};
const color=v=>v>=90?"#2358cf":v>=70?"#59c6e8":v>=50?"#f2cf5b":"#e77c47";
export function inferUniProtAccession(description="") {
  const candidate=/(?:sp|tr)\|([^|]+)\|/.exec(description)?.[1]??description.split(/\s+/)[0];
  return /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$/.test(candidate)?candidate:"";
}
export function parsePdbAtoms(text) {
  const atoms=[];let modelSeen=false;
  for(const line of String(text).split(/\r?\n/)) {
    if(line.startsWith("MODEL")){if(modelSeen)break;modelSeen=true;}
    if(line.startsWith("ENDMDL"))break;
    if(!line.startsWith("ATOM  ")&&!line.startsWith("HETATM"))continue;
    const resn=line.slice(17,20).trim();if(!AA[resn])continue;
    const alt=line[16];if(alt!==" "&&alt!=="A"&&alt!=="")continue;
    const atom=line.slice(12,16).trim();
    atoms.push({index:atoms.length,atom,resn,chain:line[21].trim(),resi:Number(line.slice(22,26)),icode:line[26].trim(),x:Number(line.slice(30,38)),y:Number(line.slice(38,46)),z:Number(line.slice(46,54)),b:parseFloat(line.slice(60,66)),elem:line.slice(76,78).trim()||atom.replace(/[0-9]/g,"")[0]});
  }return atoms;
}
export function parsePdbPlddt(text) {
  const residues=new Map();
  for(const a of parsePdbAtoms(text))if(Number.isFinite(a.b))residues.set(`${a.chain}:${a.resi}:${a.icode}`,a.b);
  const values=[...residues.values()];return values.length?{mean:values.reduce((a,b)=>a+b,0)/values.length,residueScores:residues}:null;
}
export function chainsFromAtoms(atoms) {
  const chains=new Map();
  for(const atom of atoms) {
    if(!AA[atom.resn]||[atom.x,atom.y,atom.z].some(v=>!Number.isFinite(v)))continue;
    const chain=String(atom.chain??"");if(!chains.has(chain))chains.set(chain,new Map());
    const key=`${atom.resi}:${atom.icode??""}`;const map=chains.get(chain);
    if(!map.has(key))map.set(key,{key,chain,resi:atom.resi,icode:atom.icode??"",aa:AA[atom.resn],atoms:[]});
    const r=map.get(key);if(!r.atoms.some(a=>a.atom===atom.atom))r.atoms.push(atom);
  }
  return [...chains].map(([id,residues])=>({id,residues:[...residues.values()].sort((a,b)=>a.resi-b.resi||a.icode.localeCompare(b.icode))}));
}
export function mapChain(sequence,chain,minimumResidues=30) {
  const observed=chain.residues.map(r=>r.aa).join("");
  if(observed.length<minimumResidues)return {status:"unmapped",reason:`Only ${observed.length} mapped residues; at least ${minimumResidues} are required.`};
  // A unique contiguous identity mapping also handles arbitrary author numbering/insertion codes.
  const first=sequence.indexOf(observed);
  let positions=null,method;
  if(first>=0 && sequence.indexOf(observed,first+1)<0){positions=chain.residues.map((_,i)=>first+i+1);method="unique exact sequence alignment";}
  else if(chain.residues.every(r=>!r.icode)) {
    // Missing coordinates may be retained if one numbering offset explains every observed residue.
    const reference=chain.residues[0];const candidates=[];
    for(let i=0;i<sequence.length;i++)if(sequence[i]===reference.aa){
      const offset=i+1-reference.resi;
      if(chain.residues.every(r=>sequence[r.resi+offset-1]===r.aa))candidates.push(offset);
    }
    if(candidates.length===1){positions=chain.residues.map(r=>r.resi+candidates[0]);method="100% observed identity with a unique numbering offset (gaps allowed)";}
  }
  if(!positions)return {status:"unmapped",reason:"No unique exact mapping to the submitted sequence. Check chain, accession, allele, missing segments, or sequence version."};
  const residues=chain.residues.map((r,i)=>({...r,position:positions[i]}));
  return {status:"mapped",method,chain:chain.id,residues,coverage:residues.length/sequence.length,byPosition:new Map(residues.map(r=>[r.position,r]))};
}
export function measureInterface(mapping,result,definition,{confidenceIsPlddt=false}={}) {
  if(mapping.status!=="mapped")return {status:"unavailable",reason:mapping.reason};
  const {whd,heuristicDomains,boundaryVerified,lrrEndAnnotated}=result.domainEvidence;
  if(!whd || !boundaryVerified || !lrrEndAnnotated)return {status:"unavailable",reason:"Supply source-supported WHD start/end and LRR start/end coordinates. Sequence-anchor estimates cannot define a measured WHD–LRR interface."};
  const lrr=heuristicDomains.lrr;
  const whdResidues=mapping.residues.filter(r=>r.position>=whd.start&&r.position<=whd.end);
  const lrrResidues=mapping.residues.filter(r=>r.position>=lrr.start&&r.position<=lrr.end);
  if(!whdResidues.length||!lrrResidues.length)return {status:"unavailable",reason:"The selected model has no mapped coordinates in one of the annotated domains."};
  const cutoff=definition.contactCutoffAngstrom,cut2=cutoff*cutoff;
  const heavy=r=>r.atoms.filter(a=>!['H','D'].includes(String(a.elem).toUpperCase()));
  const grid=new Map(),cell=(x,y,z)=>`${x},${y},${z}`;
  for(const r of whdResidues)for(const atom of heavy(r)){
    const key=cell(Math.floor(atom.x/cutoff),Math.floor(atom.y/cutoff),Math.floor(atom.z/cutoff));
    if(!grid.has(key))grid.set(key,[]);grid.get(key).push({atom,residue:r});
  }
  const pairs=new Map();
  for(const r of lrrResidues)for(const atom of heavy(r)){
    const cx=Math.floor(atom.x/cutoff),cy=Math.floor(atom.y/cutoff),cz=Math.floor(atom.z/cutoff);
    for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++)for(const w of grid.get(cell(cx+dx,cy+dy,cz+dz))??[]){
      if(Math.abs(r.position-w.residue.position)<definition.minimumSequenceSeparation)continue;
      const d2=(atom.x-w.atom.x)**2+(atom.y-w.atom.y)**2+(atom.z-w.atom.z)**2;
      if(d2>cut2)continue;
      const key=`${w.residue.position}:${r.position}`,prior=pairs.get(key);
      if(!prior||d2<prior.distance**2)pairs.set(key,{whdPosition:w.residue.position,whdAA:w.residue.aa,lrrPosition:r.position,lrrAA:r.aa,distance:Math.sqrt(d2),whdAuthor:`${w.residue.chain}:${w.residue.resi}${w.residue.icode}`,lrrAuthor:`${r.chain}:${r.resi}${r.icode}`,atom1:w.atom.atom,atom2:atom.atom});
    }
  }
  const contacts=[...pairs.values()].sort((a,b)=>a.lrrPosition-b.lrrPosition||a.whdPosition-b.whdPosition);
  const local=result.interfaceEvidence.localHits.map(hit=>{
    const windowPairs=contacts.filter(p=>p.lrrPosition>=hit.start&&p.lrrPosition<=hit.end);
    return {...hit,mappedResidues:mapping.residues.filter(r=>r.position>=hit.start&&r.position<=hit.end).length,windowPairs:windowPairs.length,
      constrainedPairs:windowPairs.filter(p=>hit.constrainedPositions.includes(p.lrrPosition)).length,
      pairs:windowPairs};
  });
  const confidence=residues=>{
    if(!confidenceIsPlddt)return null;
    const v=residues.map(r=>(r.atoms.find(a=>a.atom==="CA")??r.atoms[0]).b).filter(n=>Number.isFinite(n)&&n>=0&&n<=100);
    return v.length?v.reduce((a,b)=>a+b,0)/v.length:null;
  };
  return {status:"measured",totalPairs:contacts.length,contacts,localMotifs:local,cutoff,chain:mapping.chain,mappingMethod:mapping.method,
    whdCoverage:whdResidues.length/(whd.end-whd.start+1),lrrCoverage:lrrResidues.length/(lrr.end-lrr.start+1),
    whdPlddt:confidence(whdResidues),lrrPlddt:confidence(lrrResidues),
    reason:"Geometric contact counts only; no energetic stabilization or TS/TT threshold is inferred. PAE and conformational-state comparison are not available."};
}
export async function fetchAlphaFoldPdb(accession,sequence) {
  const clean=String(accession).trim().toUpperCase();
  if(!inferUniProtAccession(clean))throw new Error("Enter a valid UniProt accession.");
  const response=await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(clean)}`);
  if(!response.ok)throw new Error(`No AlphaFold DB model found for ${clean}.`);
  const data=await response.json();const prediction=data.find(p=>p.sequence===sequence)??data[0];
  if(!prediction?.pdbUrl)throw new Error("No PDB model is available.");
  if(prediction.sequence && prediction.sequence!==sequence)throw new Error("AlphaFold sequence differs from this input. Upload the appropriate allele/model instead; automatic highlighting was blocked.");
  const pdb=await fetch(prediction.pdbUrl);if(!pdb.ok)throw new Error("AlphaFold coordinates could not be retrieved.");
  return {text:await pdb.text(),format:"pdb",source:`AlphaFold DB · ${clean}`,confidenceIsPlddt:true};
}
export function createStructureViewer(element) {
  if(!globalThis.$3Dmol)throw new Error("3Dmol.js is unavailable. Reload the page.");
  const viewer=globalThis.$3Dmol.createViewer(element,{backgroundColor:"#ffffff"});let atoms=[],mapping=null;
  const select=positions=>({index:mapping?.residues.filter(r=>positions.has(r.position)).flatMap(r=>r.atoms.map(a=>a.index))??[]});
  return {
    load(text,format){viewer.clear();const model=viewer.addModel(text,format,{multimodel:false});atoms=model.selectedAtoms({});if(!atoms.length)throw new Error("No readable protein atoms in this file.");return chainsFromAtoms(atoms);},
    display(result,chain,definition,confidenceIsPlddt=false){
      mapping=mapChain(result.sequence,chain,definition.minimumMappedResidues);
      viewer.setStyle({},{cartoon:confidenceIsPlddt?{colorfunc:a=>color(a.b)}:{color:"#8aa89f"}});
      const span=(s,e)=>new Set(Array.from({length:Math.max(0,e-s+1)},(_,i)=>s+i));
      if(mapping.status==="mapped"){
        const m=result.bestMatch;
        if(m){viewer.setStyle(select(span(Math.max(1,m.start-50),Math.min(result.length,m.end+50))),{cartoon:{color:"#4bb6c6"}});
          viewer.setStyle(select(span(m.start,m.end)),{cartoon:{color:"#57e6b0"},stick:{color:"#57e6b0"}});
          viewer.setStyle(select(new Set([m.start+3])),{stick:{color:"#ffc96b"},sphere:{color:"#ffc96b",scale:.4}});}
        const positions=new Set(result.interfaceEvidence.localHits.flatMap(h=>h.constrainedPositions));
        viewer.setStyle(select(positions),{cartoon:{color:"#d68cff"},stick:{color:"#d68cff"}});
      }
      viewer.zoomTo();viewer.render();
      return {mapping,measurements:measureInterface(mapping,result,definition,{confidenceIsPlddt})};
    },
    focusResidue(position){if(mapping?.status!=="mapped")return;const sel=select(new Set([position]));viewer.addStyle(sel,{stick:{color:"#d85822"},sphere:{scale:.3,color:"#d85822"}});viewer.zoomTo(sel);viewer.render();},
    focus(result){const m=result.bestMatch;if(mapping?.status!=="mapped"||!m)return;const positions=new Set(Array.from({length:Math.min(result.length,m.end+50)-Math.max(1,m.start-50)+1},(_,i)=>Math.max(1,m.start-50)+i));viewer.zoomTo(select(positions));viewer.render();},
    reset(){viewer.zoomTo();viewer.render();},
    clear(){viewer.clear();viewer.render();mapping=null;},resize(){viewer.resize();viewer.render();}
  };
}
