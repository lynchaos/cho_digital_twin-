import { useMemo, useState } from "react";
import { runSimulation, defaultConfig } from "@/lib/simulator";
import { runFBA, computeFluxPCA, type FBAResult, type FBAFluxes } from "@/lib/fba-solver";
import { SVG_METABOLITES, type RxnId } from "@/lib/cho-network";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend,
  ResponsiveContainer, ScatterChart, Scatter, Cell,
} from "recharts";

// ── Color categories ──────────────────────────────────────────────────────────
const RXN_COLOR: Record<string, string> = {
  EX_GLC:"#2563eb", EX_LAC:"#2563eb", EX_GLN:"#2563eb", EX_GLU:"#2563eb",
  LDH_f:"#d97706",  LDH_r:"#d97706",  PDH:"#d97706",
  CS:"#16a34a",     ICDH:"#16a34a",   AKGDH:"#16a34a",
  SDH_FH:"#16a34a", MDH:"#16a34a",
  ME:"#ea580c",     PCX:"#dc2626",
  GLS:"#9333ea",    GLUD:"#9333ea",
};

// ── Arrow definitions ─────────────────────────────────────────────────────────
interface ArrowDef { id: RxnId; from: string; to: string; perpOff?: number; }
const ARROWS: ArrowDef[] = [
  { id:"EX_GLC",  from:"GLC_EXT", to:"PYR"     },
  { id:"EX_LAC",  from:"LAC",     to:"LAC_EXT"  },
  { id:"EX_GLN",  from:"GLN_EXT", to:"GLN"      },
  { id:"EX_GLU",  from:"GLU_EXT", to:"GLU"      },
  { id:"LDH_f",   from:"PYR",     to:"LAC",     perpOff: 9  },
  { id:"LDH_r",   from:"LAC",     to:"PYR",     perpOff: 9  },
  { id:"PDH",     from:"PYR",     to:"AcCoA"    },
  { id:"CS",      from:"AcCoA",   to:"CIT"      },
  { id:"ICDH",    from:"CIT",     to:"AKG"      },
  { id:"AKGDH",   from:"AKG",     to:"SUCC"     },
  { id:"SDH_FH",  from:"SUCC",    to:"MAL"      },
  { id:"MDH",     from:"MAL",     to:"OAA"      },
  { id:"ME",      from:"MAL",     to:"PYR"      },
  { id:"PCX",     from:"PYR",     to:"OAA"      },
  { id:"GLS",     from:"GLN",     to:"GLU"      },
  { id:"GLUD",    from:"GLU",     to:"AKG"      },
];

const NODE_LABELS: Record<string, string> = {
  GLC_EXT:"Glc", PYR:"PYR",  LAC:"LAC",  LAC_EXT:"Lac",
  AcCoA:"AcCoA", OAA:"OAA",  CIT:"CIT",  AKG:"α-KG",
  SUCC:"SUCC",   MAL:"MAL",  GLN_EXT:"Gln",
  GLN:"GLN",     GLU:"GLU",  GLU_EXT:"Glu",
};

const NODE_POS = Object.fromEntries(
  SVG_METABOLITES.map(m => [m.id, { x:m.x, y:m.y, ext:!!m.extracellular }])
);

// ── SVG geometry helpers ──────────────────────────────────────────────────────
function arrowGeom(
  from:{x:number;y:number}, to:{x:number;y:number},
  fromR:number, toR:number, perpOff:number
) {
  const dx=to.x-from.x, dy=to.y-from.y;
  const len=Math.sqrt(dx*dx+dy*dy);
  if(len<1) return null;
  const nx=dx/len, ny=dy/len, px=-ny, py=nx;
  const x1=from.x+px*perpOff+nx*fromR,   y1=from.y+py*perpOff+ny*fromR;
  const tipX=to.x+px*perpOff-nx*toR,     tipY=to.y+py*perpOff-ny*toR;
  const x2=tipX-nx*7,                     y2=tipY-ny*7;
  const ax=tipX-nx*9-py*4.5,             ay=tipY-ny*9+px*4.5;
  const bx=tipX-nx*9+py*4.5,             by=tipY-ny*9-px*4.5;
  return {x1,y1,x2,y2,tipX,tipY,ax,ay,bx,by};
}

function fluxW(v:number, maxF:number) {
  const f=Math.abs(v);
  return f<0.001 ? 1 : Math.max(1.5, Math.min(9, 1.5+7.5*f/Math.max(maxF,0.01)));
}

// ── FluxMap SVG component ─────────────────────────────────────────────────────
function FluxMap({ fluxes, maxFlux }: { fluxes: FBAFluxes; maxFlux: number }) {
  const [hov, setHov] = useState<string|null>(null);

  return (
    <svg viewBox="0 0 830 435" style={{width:"100%",height:"auto",display:"block"}}>

      {/* ─ Colour legend ─ */}
      {([["#2563eb","Exchange"],["#d97706","Glycolysis"],["#16a34a","TCA cycle"],
         ["#ea580c","Malic enz."],["#9333ea","Amino acids"]] as [string,string][]).map(([c,lab],i)=>(
        <g key={lab} transform={`translate(${8+i*122},420)`}>
          <rect x={0} y={-4} width={14} height={5} rx={2} fill={c}/>
          <text x={17} y={1} fontSize={9} fill="#64748b">{lab}</text>
        </g>
      ))}

      {/* ─ Reaction arrows ─ */}
      {ARROWS.map(arw=>{
        const fn=NODE_POS[arw.from], tn=NODE_POS[arw.to];
        if(!fn||!tn) return null;
        const flux=fluxes[arw.id]??0;
        const color=RXN_COLOR[arw.id]??"#94a3b8";
        const w=fluxW(flux,maxFlux);
        const isHov=hov===arw.id;
        const alpha=flux<0.001 ? 0.13 : (isHov ? 1 : 0.72);
        const fromR=fn.ext?20:26, toR=tn.ext?20:26;
        const geo=arrowGeom(fn,tn,fromR,toR,arw.perpOff??0);
        if(!geo) return null;

        return (
          <g key={arw.id} opacity={alpha} style={{cursor:"pointer"}}
             onMouseEnter={()=>setHov(arw.id)} onMouseLeave={()=>setHov(null)}>
            <line x1={geo.x1} y1={geo.y1} x2={geo.x2} y2={geo.y2}
              stroke={color} strokeWidth={isHov?w+2:w} strokeLinecap="round"/>
            <path d={`M${geo.ax},${geo.ay}L${geo.tipX},${geo.tipY}L${geo.bx},${geo.by}`}
              stroke={color} strokeWidth={isHov?2.6:1.8}
              fill="none" strokeLinecap="round" strokeLinejoin="round"/>
            {isHov && (
              <g>
                <rect x={(geo.x1+geo.x2)/2-26} y={(geo.y1+geo.y2)/2-11}
                  width={52} height={16} rx={3}
                  fill="white" stroke={color} strokeWidth={1.2} opacity={0.96}/>
                <text x={(geo.x1+geo.x2)/2} y={(geo.y1+geo.y2)/2+2}
                  textAnchor="middle" fontSize={9.5} fill={color} fontWeight="600">
                  {flux.toFixed(4)}
                </text>
              </g>
            )}
          </g>
        );
      })}

      {/* ─ Metabolite nodes ─ */}
      {SVG_METABOLITES.map(met=>{
        const isExt=!!met.extracellular;
        const r=isExt?20:26;
        const label=NODE_LABELS[met.id]??met.id;
        return (
          <g key={met.id} transform={`translate(${met.x},${met.y})`}>
            <circle r={r}
              fill={isExt?"#eff6ff":"white"}
              stroke={isExt?"#93c5fd":"#334155"}
              strokeWidth={isExt?1.5:2.2}
              strokeDasharray={isExt?"4 3":undefined}/>
            <text textAnchor="middle" dominantBaseline="middle"
              fontSize={isExt?8.5:9.5} fontWeight={isExt?"400":"600"}
              fill={isExt?"#1e40af":"#1e293b"}>
              {label}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

// ── Flux info groups for side-panel table ─────────────────────────────────────
const FLUX_GROUPS: { title:string; rows:{id:RxnId;label:string;desc:string}[] }[] = [
  { title:"Glycolysis / Pyruvate", rows:[
    {id:"EX_GLC" as RxnId, label:"EX_GLC", desc:"Glc → 2 Pyr (lump)"},
    {id:"PDH"    as RxnId, label:"PDH",    desc:"Pyr → AcCoA + CO₂"},
    {id:"LDH_f"  as RxnId, label:"LDH →",  desc:"Pyr → Lac"},
    {id:"LDH_r"  as RxnId, label:"LDH ←",  desc:"Lac → Pyr"},
    {id:"EX_LAC" as RxnId, label:"EX_LAC", desc:"Net lactate exchange"},
  ]},
  { title:"TCA Cycle", rows:[
    {id:"CS"     as RxnId, label:"CS",     desc:"OAA+AcCoA → Cit"},
    {id:"ICDH"   as RxnId, label:"ICDH",   desc:"Cit → α-KG + CO₂"},
    {id:"AKGDH"  as RxnId, label:"αKGDH",  desc:"α-KG → Succ + CO₂"},
    {id:"SDH_FH" as RxnId, label:"SDH/FH", desc:"Succ → Mal"},
    {id:"MDH"    as RxnId, label:"MDH",    desc:"Mal → OAA"},
  ]},
  { title:"Anaplerosis", rows:[
    {id:"ME"  as RxnId, label:"ME",  desc:"Mal → Pyr + CO₂"},
    {id:"PCX" as RxnId, label:"PCX", desc:"Pyr + CO₂ → OAA"},
  ]},
  { title:"Amino Acids", rows:[
    {id:"EX_GLN" as RxnId, label:"EX_GLN", desc:"Gln uptake from broth"},
    {id:"GLS"    as RxnId, label:"GLS",    desc:"Gln → Glu + NH₃"},
    {id:"EX_GLU" as RxnId, label:"EX_GLU", desc:"Glu uptake from broth"},
    {id:"GLUD"   as RxnId, label:"GLUD",   desc:"Glu → α-KG + NH₃"},
  ]},
];

// ── Utility ───────────────────────────────────────────────────────────────────
function dayColor(day:number, maxDay:number):string {
  const t=Math.max(0,Math.min(1,day/maxDay));
  return `rgb(${Math.round(59+t*161)},${Math.round(130*(1-t))},${Math.round(246-t*196)})`;
}

// ── Main page ─────────────────────────────────────────────────────────────────
export default function PcdFBAPage() {
  const simResults = useMemo(()=>runSimulation(defaultConfig()),[]);

  const fbaResults = useMemo(()=>
    simResults.map(tp=>runFBA(tp.q_Glc, tp.q_Lac, tp.q_Gln, tp.q_Glu)),
  [simResults]);

  const pca = useMemo(()=>computeFluxPCA(fbaResults),[fbaResults]);

  const [tIdx, setTIdx] = useState<number>(Math.floor(simResults.length/2));

  const cur:FBAResult|undefined = fbaResults[tIdx];
  const curT:number = simResults[tIdx]?.t ?? 0;
  const maxDay:number = simResults[simResults.length-1]?.t ?? 14;

  const maxFlux = useMemo(()=>{
    let m=0.1;
    for(const r of fbaResults)
      for(const v of Object.values(r.fluxes))
        if(Math.abs(v as number)>m) m=Math.abs(v as number);
    return m;
  },[fbaResults]);

  type TSRow = {t:number;EX_GLC:number;PDH:number;CS:number;AKGDH:number;
                GLS:number;GLUD:number;LDH_f:number;EX_LAC:number};
  const tsData = useMemo(():TSRow[]=>{
    return simResults.map((tp,i)=>{
      const f=fbaResults[i]?.fluxes;
      if(!f) return null as unknown as TSRow;
      return {t:+tp.t.toFixed(2), EX_GLC:f.EX_GLC, PDH:f.PDH, CS:f.CS,
              AKGDH:f.AKGDH, GLS:f.GLS, GLUD:f.GLUD, LDH_f:f.LDH_f, EX_LAC:f.EX_LAC};
    }).filter(Boolean);
  },[simResults,fbaResults]);

  const pcaData = useMemo(()=>pca.scores.map(([pc1,pc2],i)=>({
    pc1:+pc1.toFixed(5), pc2:+pc2.toFixed(5), day:+(simResults[i]?.t??0).toFixed(1),
  })),[pca,simResults]);

  const loadingData = useMemo(()=>pca.rxnOrder.map((rxn,i)=>({
    rxn, pc1:+(pca.loadings[0][i]??0).toFixed(4), pc2:+(pca.loadings[1][i]??0).toFixed(4),
  })),[pca]);

  return (
    <div className="pcdfba-page">

      {/* ── Header ─────────────────────────────────────────────────────── */}
      <div className="pcdfba-header">
        <div className="pcdfba-header-left">
          <h1 className="pcdfba-title">PC-dFBA Intracellular Flux Analysis</h1>
          <p className="pcdfba-subtitle">
            Condensed 10-metabolite / 16-reaction CHO central carbon metabolism
            model (glycolysis → TCA → amino acid catabolism).
            Exchange rates from the ODE constrain FBA; most fluxes are uniquely
            determined by mass-balance + parsimony.
          </p>
        </div>
        <div className="pcdfba-alg-box">
          <div className="pcdfba-alg-title">Analytical FBA solution (Eqs. 27–33 condensed)</div>
          <code className="pcdfba-alg-code">
            GLUD = EX_GLN + EX_GLU &nbsp;&nbsp; PDH = 2·EX_GLC − EX_LAC + GLUD<br/>
            CS = ICDH = PDH &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; AKGDH = SDH = PDH + GLUD<br/>
            ME = GLUD (parsimony) &nbsp;&nbsp; MDH = PDH &nbsp;&nbsp; PCX = 0
          </code>
        </div>
      </div>

      {/* ── Main: flux map + info panel ────────────────────────────────── */}
      <div className="pcdfba-main">

        {/* Flux map */}
        <div className="pcdfba-svg-col">
          <div className="pcdfba-card">
            <div className="pcdfba-card-title">
              Metabolic Flux Map
              <span className="pcdfba-day-badge">Day {curT.toFixed(1)}</span>
            </div>
            {cur && <FluxMap fluxes={cur.fluxes} maxFlux={maxFlux}/>}
            <div className="pcdfba-slider-row">
              <span className="pcdfba-slider-label">Day 0</span>
              <input type="range" className="pcdfba-slider"
                min={0} max={simResults.length-1} value={tIdx}
                onChange={e=>setTIdx(+e.target.value)}/>
              <span className="pcdfba-slider-label">Day {maxDay.toFixed(0)}</span>
            </div>
          </div>
        </div>

        {/* Info panel */}
        <div className="pcdfba-info-col">
          {cur && FLUX_GROUPS.map(grp=>(
            <div key={grp.title} className="pcdfba-card pcdfba-flux-group">
              <div className="pcdfba-flux-group-title">{grp.title}</div>
              <table className="pcdfba-flux-table">
                <tbody>
                  {grp.rows.map(row=>{
                    const v=cur.fluxes[row.id]??0;
                    const col=RXN_COLOR[row.id]??"#475569";
                    return (
                      <tr key={row.id}>
                        <td className="pcdfba-rxn-label" style={{color:col}}>{row.label}</td>
                        <td className="pcdfba-rxn-val">{v.toFixed(4)}</td>
                        <td className="pcdfba-rxn-unit">mmol/Mc/d</td>
                        <td className="pcdfba-rxn-desc">{row.desc}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          ))}

          {cur && (
            <div className="pcdfba-card">
              <span className={`pcdfba-status-badge ${cur.status==="feasible"
                ?"pcdfba-status-ok":"pcdfba-status-warn"}`}>
                {cur.status==="feasible" ? "✓ Feasible" : "⚠ Infeasible"}
              </span>
              {cur.warnings.length>0 && (
                <ul className="pcdfba-warnings">
                  {cur.warnings.map((w,i)=><li key={i}>{w}</li>)}
                </ul>
              )}
            </div>
          )}
        </div>
      </div>

      {/* ── Flux time series ───────────────────────────────────────────── */}
      <div className="pcdfba-charts-row">
        <div className="pcdfba-card pcdfba-chart-card">
          <div className="pcdfba-card-title">Glycolysis / TCA Fluxes over Culture</div>
          <ResponsiveContainer width="100%" height={200}>
            <LineChart data={tsData} margin={{top:4,right:10,left:-10,bottom:0}}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0"/>
              <XAxis dataKey="t" tick={{fontSize:10}}
                label={{value:"Day",position:"insideBottomRight",offset:-4,fontSize:10}}/>
              <YAxis tick={{fontSize:10}}/>
              <Tooltip formatter={(v:unknown)=>(+(v as number)).toFixed(4)}
                       labelFormatter={l=>`Day ${l}`}/>
              <Legend wrapperStyle={{fontSize:"10px"}}/>
              <Line type="monotone" dataKey="EX_GLC" name="EX_GLC"  stroke="#2563eb" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="PDH"    name="PDH"     stroke="#d97706" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="CS"     name="CS"      stroke="#16a34a" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="AKGDH"  name="αKGDH"   stroke="#0891b2" dot={false} strokeWidth={1.8}/>
            </LineChart>
          </ResponsiveContainer>
        </div>

        <div className="pcdfba-card pcdfba-chart-card">
          <div className="pcdfba-card-title">Amino Acid / Lactate Fluxes over Culture</div>
          <ResponsiveContainer width="100%" height={200}>
            <LineChart data={tsData} margin={{top:4,right:10,left:-10,bottom:0}}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0"/>
              <XAxis dataKey="t" tick={{fontSize:10}}
                label={{value:"Day",position:"insideBottomRight",offset:-4,fontSize:10}}/>
              <YAxis tick={{fontSize:10}}/>
              <Tooltip formatter={(v:unknown)=>(+(v as number)).toFixed(4)}
                       labelFormatter={l=>`Day ${l}`}/>
              <Legend wrapperStyle={{fontSize:"10px"}}/>
              <Line type="monotone" dataKey="GLS"    name="GLS"    stroke="#9333ea" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="GLUD"   name="GLUD"   stroke="#a855f7" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="LDH_f"  name="LDH→"   stroke="#d97706" dot={false} strokeWidth={1.8}/>
              <Line type="monotone" dataKey="EX_LAC" name="EX_LAC" stroke="#2563eb" dot={false} strokeWidth={1.8}/>
            </LineChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* ── PC Analysis ────────────────────────────────────────────────── */}
      <div className="pcdfba-charts-row pcdfba-pca-row">
        <div className="pcdfba-card pcdfba-chart-card pcdfba-pca-card">
          <div className="pcdfba-card-title">
            Metabolic State Trajectory (PC Space)
            <span className="pcdfba-small-note">
              PC1: {(pca.explained[0]*100).toFixed(0)}% var &nbsp;·&nbsp;
              PC2: {(pca.explained[1]*100).toFixed(0)}% var
            </span>
          </div>
          <div className="pcdfba-pca-legend">
            <div className="pcdfba-pca-grad"/>
            <span className="pcdfba-pca-leg-label">Day 0</span>
            <span className="pcdfba-pca-leg-label" style={{marginLeft:"auto"}}>Day {maxDay.toFixed(0)}</span>
          </div>
          {pcaData.length>3 ? (
            <ResponsiveContainer width="100%" height={210}>
              <ScatterChart margin={{top:5,right:15,left:-15,bottom:5}}>
                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0"/>
                <XAxis type="number" dataKey="pc1" name="PC1"
                  tick={{fontSize:10}} label={{value:"PC1",position:"insideBottom",offset:-2,fontSize:10}}/>
                <YAxis type="number" dataKey="pc2" name="PC2"
                  tick={{fontSize:10}} label={{value:"PC2",angle:-90,position:"insideLeft",fontSize:10}}/>
                <Tooltip cursor={{strokeDasharray:"3 3"}}
                  formatter={(v:unknown,n:unknown)=>[(+(v as number)).toFixed(5),n as string]}
                  labelFormatter={(_,pl)=>`Day ${(pl as {payload:{day:number}}[])?.[0]?.payload?.day??""}`}/>
                <Scatter data={pcaData} name="State">
                  {pcaData.map((d,i)=><Cell key={i} fill={dayColor(d.day,maxDay)}/>)}
                </Scatter>
              </ScatterChart>
            </ResponsiveContainer>
          ):(
            <div className="pcdfba-pca-empty">Insufficient data for PCA</div>
          )}
        </div>

        <div className="pcdfba-card pcdfba-chart-card">
          <div className="pcdfba-card-title">PC Loadings — Reaction Contributions</div>
          <p className="pcdfba-chart-sub">
            PC1 (blue) and PC2 (orange) loading vectors: how each reaction drives
            metabolic-state variation across the culture trajectory.
          </p>
          <ResponsiveContainer width="100%" height={220}>
            <LineChart data={loadingData} margin={{top:4,right:10,left:-10,bottom:36}}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0"/>
              <XAxis dataKey="rxn" tick={{fontSize:8}} angle={-40} textAnchor="end" interval={0}/>
              <YAxis tick={{fontSize:9}}/>
              <Tooltip/>
              <Legend wrapperStyle={{fontSize:"10px"}}/>
              <Line type="monotone" dataKey="pc1" name="PC1" stroke="#2563eb" dot={{r:3}} strokeWidth={1.5}/>
              <Line type="monotone" dataKey="pc2" name="PC2" stroke="#f59e0b" dot={{r:3}} strokeWidth={1.5}/>
            </LineChart>
          </ResponsiveContainer>
        </div>
      </div>

    </div>
  );
}
