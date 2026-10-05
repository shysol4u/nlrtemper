import {SHORT_STRETCHES} from './short-reference.mjs';
// Shared columns transcribed from the user's supplied manuscript figure.
// Shared gap columns are retained, never counted as amino-acid coordinates.
export const REFERENCE_ROWS=['NLEELDLVGC-KSLV-TLPSSIQNA','NLEYVNLYQC-SNLE-EVHHSLGCC'];
export function jointReferenceComparison(record,references){
 const refs=['SNC1','N'].map((id,k)=>{const r=references.find(r=>r.id===id),start=r?.sequence.indexOf(SHORT_STRETCHES[id])??-1;return {id,start:start+1,end:start+SHORT_STRETCHES[id].length,row:REFERENCE_ROWS[k]};});
 if(refs.some(r=>!r.start))return {error:'Both reference stretches are required for the joint alignment.'};
 const query=record.sequence,n=REFERENCE_ROWS[0].length,m=query.length,w=m+1,s=new Int32Array((n+1)*w),t=new Uint8Array(s.length);
 const column=i=>REFERENCE_ROWS.map(r=>r[i]),skip=i=>column(i).filter(c=>c!=='-').length*-2;
 for(let i=1;i<=n;i++){s[i*w]=s[(i-1)*w]+skip(i-1);t[i*w]=1;}
 for(let i=1;i<=n;i++)for(let j=1;j<=m;j++){
  const k=i*w+j,cs=column(i-1),d=s[k-w-1]+cs.reduce((v,c)=>v+(c==='-'?-2:c===query[j-1]?2:-1),0),u=s[k-w]+skip(i-1),l=s[k-1]-4;
  s[k]=Math.max(d,u,l);t[k]=cs.every(c=>c==='-')&&s[k]===u?1:s[k]===d?0:s[k]===u?1:2;
 }
 let end=0,score=s[n*w],ties=0;for(let j=1;j<=m;j++){if(s[n*w+j]>score){end=j;score=s[n*w+j];ties=0;}else if(s[n*w+j]===score)ties++;}
 let i=n,j=end;const columns=[];
 while(i>0){const d=t[i*w+j];if(d===0&&j){columns.push({referenceIndex:i-1,refs:column(--i),query:query[--j],queryPosition:j+1});}else if(d===1||!j){columns.push({referenceIndex:i-1,refs:column(--i),query:'-',queryPosition:null});}else{columns.push({referenceIndex:null,refs:['-','-'],query:query[--j],queryPosition:j+1});}}
 columns.reverse();const positions=refs.map(r=>r.start-1);
 for(const c of columns){c.refPositions=c.refs.map((aa,k)=>aa==='-'?null:++positions[k]);c.motif=c.referenceIndex!=null&&c.referenceIndex<5;c.second=c.referenceIndex===3;}
 const stats=refs.map((r,k)=>{const paired=columns.filter(c=>c.refs[k]!=='-'&&c.query!=='-').length,identical=columns.filter(c=>c.refs[k]!=='-'&&c.refs[k]===c.query).length;return {identity:identical/23,coverage:paired/23};});
 return {refs,columns,start:j+1,end,score,tied:ties>0,stats,mapped:columns.find(c=>c.second)?.queryPosition??null,quality:stats.some(s=>s.identity>=.5&&s.coverage>=.8)?'Exploratory joint correspondence':'Low sequence agreement'};
}
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export function jointReferenceHtml(record){
 const c=record.jointComparison;if(!c)return '';if(c.error)return `<p>${esc(c.error)}</p>`;
 const row=(label,start,end,k)=>`<div class="short-line joint-line"><span>${label} ${start} </span><code>${c.columns.map(x=>{const aa=k===2?x.query:x.refs[k],p=k===2?x.queryPosition:x.refPositions[k];return `<span class="${x.second?'second-e':x.motif?'motif':''}" title="${esc(p==null?'Gap':`${label} ${aa}${p}`)}">${aa}</span>`;}).join('')}</code><span> ${end}</span></div>`;
 return `<section class="short-pair joint-pair"><h4>Combined SNC1 + tobacco N + query</h4><div class="alignment" aria-label="Joint three-sequence alignment">${row('SNC1',c.refs[0].start,c.refs[0].end,0)}${row('Tobacco N',c.refs[1].start,c.refs[1].end,1)}${row('Query',c.start,c.end,2)}<div class="short-line joint-line"><span>Agreement </span><code>${c.columns.map(x=>x.query!=='-'&&x.refs.every(r=>r===x.query)?'*':x.query!=='-'&&x.refs.includes(x.query)?'|':' ').join('')}</code></div></div><p>${c.stats.map((s,k)=>`${k?'Tobacco N':'SNC1'}: ${Math.round(s.identity*100)}% identity, ${Math.round(s.coverage*100)}% paired coverage`).join(' · ')}.</p><p>E640 / Y646-equivalent in query: <strong>${c.mapped==null?'gap':esc(record.sequence[c.mapped-1])+c.mapped}</strong>. ${esc(c.quality)}${c.tied?' · equally scoring placements exist':''}.</p><p class="tiny">SNC1 and tobacco N retain the shared columns in your supplied figure. The query is fitted against both references with equal weight; its position can differ from either individual comparison. * = all three agree; | = query agrees with one reference. Gaps do not consume residue numbers. This is sequence correspondence, not structural superposition or a temperature prediction.</p></section>`;
}
