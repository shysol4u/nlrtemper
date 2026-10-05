// Window alignment is exploratory correspondence, never a phenotype transfer.
export function alignWindows(a,b){
 const n=a.length,m=b.length,w=m+1,s=new Int32Array((n+1)*w),t=new Uint8Array(s.length);
 for(let i=1;i<=n;i++){s[i*w]=-2*i;t[i*w]=1;}for(let j=1;j<=m;j++){s[j]=-2*j;t[j]=2;}
 for(let i=1;i<=n;i++)for(let j=1;j<=m;j++){const k=i*w+j,v=[s[k-w-1]+(a[i-1]===b[j-1]?2:-1),s[k-w]-2,s[k-1]-2];s[k]=Math.max(...v);t[k]=v.indexOf(s[k]);}
 let i=n,j=m,x='',y='';while(i||j){const d=t[i*w+j];if(d===0){x=a[--i]+x;y=b[--j]+y;}else if(d===1){x=a[--i]+x;y='-'+y;}else{x='-'+x;y=b[--j]+y;}}
 const pairs=[...x].filter((c,k)=>c!=='-'&&y[k]!=='-').length,identical=[...x].filter((c,k)=>c!=='-'&&c===y[k]).length;
 return {reference:x,query:y,score:s[n*w+m],identity:pairs?identical/pairs:0,coverage:pairs/Math.max(n,m)};
}
export function referenceComparisons(result,references){return references.map(ref=>{
 const start=Math.max(1,ref.anchor-50),end=Math.min(ref.sequence.length,ref.anchor+54),seq=ref.sequence.slice(start-1,end);
 const candidates=result.matches.map(hit=>{const queryStart=Math.max(1,hit.start-50),queryEnd=Math.min(result.length,hit.end+50),alignment=alignWindows(seq,result.sequence.slice(queryStart-1,queryEnd));let rp=start-1,qp=queryStart-1,mapped=null;for(let k=0;k<alignment.reference.length;k++){if(alignment.reference[k]!=='-')rp++;if(alignment.query[k]!=='-')qp++;if(rp===ref.anchor+3&&alignment.reference[k]!=='-'){mapped=alignment.query[k]==='-'?null:qp;break;}}return {hit,queryStart,queryEnd,alignment,mapped};}).sort((a,b)=>b.alignment.score-a.alignment.score||a.hit.start-b.hit.start);
 return {ref,start,end,best:candidates[0]??null,tied:candidates.length>1&&candidates[0].alignment.score===candidates[1].alignment.score};
});}
