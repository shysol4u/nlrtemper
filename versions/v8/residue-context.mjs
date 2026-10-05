// Approximate solvent-accessible area by sampled rolling-probe spheres, selected chain only.
const RADII={C:1.7,N:1.55,O:1.52,S:1.8,P:1.8};
const EXPECTED={A:5,R:11,N:8,D:8,C:6,Q:9,E:9,G:4,H:10,I:8,L:8,K:9,M:8,F:11,P:7,S:6,T:7,W:14,Y:12,V:7};
export function residueContext(mapping,result,{confidenceIsPlddt=false,samples=240,measurements=null}={}){
 if(mapping.status!=='mapped')return {status:'unavailable',reason:mapping.reason,rows:[]};
 const atoms=mapping.residues.flatMap(r=>r.atoms).filter(a=>!['H','D'].includes(String(a.elem).toUpperCase()));
 const radius=a=>(RADII[String(a.elem).toUpperCase()]??1.7)+1.4;
 const grid=new Map(),key=(x,y,z)=>`${x},${y},${z}`,cell=8;
 for(const a of atoms){const k=key(Math.floor(a.x/cell),Math.floor(a.y/cell),Math.floor(a.z/cell));if(!grid.has(k))grid.set(k,[]);grid.get(k).push(a);}
 const points=Array.from({length:samples},(_,i)=>{const y=1-2*(i+.5)/samples,r=Math.sqrt(1-y*y),phi=i*Math.PI*(3-Math.sqrt(5));return [Math.cos(phi)*r,y,Math.sin(phi)*r];});
 const positions=[...new Set(result.matches.flatMap(m=>Array.from({length:m.end-m.start+1},(_,i)=>m.start+i)))].sort((a,b)=>a-b);
 const rows=positions.map(position=>{const r=mapping.byPosition.get(position),motifs=result.matches.filter(m=>position>=m.start&&position<=m.end).map(m=>`${m.sequence}:${m.start}-${m.end}`);if(!r)return {position,aa:result.sequence[position-1],motifs,secondary:'Unmapped',sasa:null,confidence:null,contacts:null};
  let area=0;for(const a of r.atoms.filter(a=>!['H','D'].includes(String(a.elem).toUpperCase()))){const ra=radius(a),cx=Math.floor(a.x/cell),cy=Math.floor(a.y/cell),cz=Math.floor(a.z/cell),neighbors=[];
   for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++)for(const b of grid.get(key(cx+x,cy+y,cz+z))??[]){if(a===b)continue;if((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2<(ra+radius(b))**2)neighbors.push(b);}
   let visible=0;for(const p of points){const x=a.x+p[0]*ra,y=a.y+p[1]*ra,z=a.z+p[2]*ra;if(!neighbors.some(b=>(x-b.x)**2+(y-b.y)**2+(z-b.z)**2<radius(b)**2))visible++;}area+=4*Math.PI*ra*ra*visible/samples;
  }
  const ca=r.atoms.find(a=>a.atom==='CA'),ss=ca?.ss;
  const backbone=['N','CA','C','O'].every(n=>r.atoms.some(a=>a.atom===n));
  const secondary=!backbone?'Unknown (incomplete backbone)':ss==='h'?'Helix':ss==='s'?'Strand':ss==='c'?'Coil / turn':'Unknown';
  const count=r.atoms.filter(a=>!['H','D'].includes(String(a.elem).toUpperCase())).length;
  return {position,aa:r.aa,motifs,secondary,outsideCoil:secondary==='Helix'||secondary==='Strand'?'Yes':secondary==='Coil / turn'?'No':'Unknown',sasa:area,atomsComplete:count>=(EXPECTED[r.aa]??Infinity),confidence:confidenceIsPlddt&&Number.isFinite(ca?.b)?ca.b:null,contacts:measurements?.status==='measured'?measurements.contacts.filter(c=>c.lrrPosition===position||c.whdPosition===position).length:null};
 });return {status:'measured',rows,method:'3Dmol secondary structure; approximate rolling-probe SASA (1.4 Å probe, 240 sphere samples/atom), selected chain only. Coil/turn is not proof of a flexible loop. SASA depends on missing atoms, assembly and model quality; no exposure cutoff or phenotype effect is inferred.'};
}
