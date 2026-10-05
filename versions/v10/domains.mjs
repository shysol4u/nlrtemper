// Import existing annotations; this browser does not run HMMER or InterProScan.
export function parseDomains(text,format,records){
 const output=new Map(),byId=new Map(records.map(r=>[r.id,r]));
 const add=(id,length,d)=>{const r=byId.get(id);if(!r)throw new Error(`Annotation protein ${id} does not match an analyzed FASTA identifier.`);if(Number(length)!==r.sequence.length)throw new Error(`${id}: annotation length differs from FASTA.`);if(!Number.isInteger(d.start)||!Number.isInteger(d.end)||d.start<1||d.end<d.start||d.end>r.sequence.length)throw new Error(`${id}: invalid domain coordinates.`);const a=output.get(id)??{domains:[],source:'Imported annotations (identifier and length matched; sequence identity not established)'};a.domains.push(d);output.set(id,a);};
 if(format==='json'){
  const list=JSON.parse(text);if(!Array.isArray(list))throw new Error('JSON must be an array of {id, sequence, annotation} records.');
  for(const x of list){const r=byId.get(x.id);if(!r||x.sequence!==r.sequence)throw new Error(`${x.id}: JSON must contain the exact analyzed sequence.`);if(!x.annotation?.source)throw new Error('JSON annotation requires a source.');output.set(x.id,x.annotation);}return output;
 }
 for(const line of text.split(/\r?\n/)){if(!line.trim()||line.startsWith('#'))continue;
  if(format==='interpro'){const c=line.split('\t');if(c.length<11)throw new Error('Expected InterProScan TSV with at least 11 columns.');add(c[0],c[2],{name:c[5]||c[4],accession:c[4],source:c[3]+' / InterProScan',start:Number(c[6]),end:Number(c[7]),score:c[8]});}
  else {const c=line.trim().split(/\s+/);if(c.length<22)throw new Error('Expected HMMER domtblout.');if(!Number.isFinite(Number(c[12])))throw new Error('Invalid domain E-value.');if(Number(c[12])>1e-5)continue;const scan=format==='hmmscan';add(c[scan?3:0],c[scan?5:2],{name:c[scan?0:3],accession:c[scan?1:4],start:Number(c[17]),end:Number(c[18]),source:`HMMER ${format}; independent E ≤ 1e-5`,score:c[12]});}
 }
 if(!output.size)throw new Error('No accepted annotation rows.');return output;
}
