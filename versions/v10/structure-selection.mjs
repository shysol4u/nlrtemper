// Inspection state is separate from the automatic motif ranking and phenotype call.
export function resolveStructureSelection(result,choice='preferred',start=null,end=null){
 let ranges;
 if(choice==='custom'){
  start=Number(start);end=Number(end);
  if(!Number.isInteger(start)||!Number.isInteger(end)||start<1||end<start||end>result.length)throw new Error(`Choose a residue range from 1 to ${result.length}, with end ≥ start.`);
  ranges=[{start,end,sequence:result.sequence.slice(start-1,end),custom:true}];
 }else if(choice==='all')ranges=result.matches;
 else if(choice==='preferred')ranges=result.bestMatch?[result.bestMatch]:[];
 else {const hit=result.matches.find(m=>m.start===Number(choice));if(!hit)throw new Error('This motif is no longer available. Choose another motif.');ranges=[hit];}
 const positions=new Set(),context=new Set(),seconds=new Set();
 for(const r of ranges){for(let p=r.start;p<=r.end;p++)positions.add(p);for(let p=Math.max(1,r.start-50);p<=Math.min(result.length,r.end+50);p++)context.add(p);if(!r.custom)seconds.add(r.start+3);}
 return {choice,ranges,positions,context,seconds,label:choice==='all'?`All ${ranges.length} detected motifs`:choice==='custom'?`Custom region ${start}–${end}`:ranges.length?`${ranges[0].sequence} · ${ranges[0].start}–${ranges[0].end}`:'No qualifying motif'};
}
export const authorCoordinate=r=>`${r.chain||'(blank)'}:${r.resi}${r.icode??''}`;
const residueKey=r=>`${r.chain}:${r.resi}:${r.icode??''}`;
const heavy=r=>r.atoms.filter(a=>!['H','D'].includes(String(a.elem).toUpperCase()));
export function selectionContacts(mapping,chains,selection,cutoff=4.5){
 if(mapping?.status!=='mapped')return {status:'unavailable',reason:mapping?.reason??'Load a matching structure.',contacts:[],neighbors:[]};
 const selected=mapping.residues.filter(r=>selection.positions.has(r.position)),selectedKeys=new Set(selected.map(residueKey)),mapped=new Map(mapping.residues.map(r=>[residueKey(r),r]));
 const grid=new Map(),cell=(x,y,z)=>`${x},${y},${z}`;
 for(const chain of chains)for(const raw of chain.residues){if(selectedKeys.has(residueKey(raw)))continue;const r=mapped.get(residueKey(raw))??raw;
  for(const atom of heavy(r)){const key=cell(Math.floor(atom.x/cutoff),Math.floor(atom.y/cutoff),Math.floor(atom.z/cutoff));if(!grid.has(key))grid.set(key,[]);grid.get(key).push({r,atom});}
 }
 const pairs=new Map();
 for(const r of selected)for(const a of heavy(r)){const x=Math.floor(a.x/cutoff),y=Math.floor(a.y/cutoff),z=Math.floor(a.z/cutoff);
  for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++)for(const {r:b,atom} of grid.get(cell(x+dx,y+dy,z+dz))??[]){
   const d2=(a.x-atom.x)**2+(a.y-atom.y)**2+(a.z-atom.z)**2;if(d2>cutoff**2)continue;
   const key=`${r.position}|${residueKey(b)}`,old=pairs.get(key);if(!old||d2<old.distance**2)pairs.set(key,{selected:r,neighbor:b,distance:Math.sqrt(d2),separation:b.position==null?null:Math.abs(r.position-b.position)});
  }
 }
 const contacts=[...pairs.values()].sort((a,b)=>a.selected.position-b.selected.position||a.distance-b.distance);
 return {status:'measured',contacts,neighbors:[...new Map(contacts.map(p=>[residueKey(p.neighbor),p.neighbor])).values()],cutoff,mapped:selected.length,total:selection.positions.size};
}
export function labelResidues(mapping,selection,nearby,mode='off'){
 if(mapping?.status!=='mapped'||mode==='off')return [];
 const selected=mapping.residues.filter(r=>selection.positions.has(r.position));
 const rows=mode==='all'?mapping.residues:mode==='neighbors'?[...selected,...nearby.neighbors]:selected;
 return [...new Map(rows.map(r=>[residueKey(r),r])).values()].map(r=>({residue:r,text:r.position==null?`${r.aa} [${authorCoordinate(r)}]`:`${r.aa}${r.position} [${authorCoordinate(r)}]`}));
}
