/**
 * GEM Reduction Pipeline Page
 *
 * Interactive visualization of the iCHO1766 genome-scale model reduction
 * pipeline from Antonakoudis & Richelle (2026) npj Systems Biology and
 * Applications. https://doi.org/10.1038/s41540-026-00704-4
 *
 * Live backend: FastAPI + COBRApy running iCHOv1 (6,663 reactions).
 * Uses async job-poll pattern to avoid proxy timeouts.
 */

import React, { useState, useCallback, useRef, useEffect } from "react";
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, Cell,
} from "recharts";

// ── Types ─────────────────────────────────────────────────────────────────────

interface ModelInfo {
  model_id: string;
  reactions: number;
  metabolites: number;
  genes: number;
  exchange_count: number;
  sample_exchanges: string[];
  core_exchanges_present: string[];
}

interface StepResult {
  step: number;
  name: string;
  reactions: number;
  metabolites: number;
  genes: number;
  exchanges_kept?: number;
  demand_reactions?: number;
  infeasible_timepoints?: number;
  duration_s?: number;
  note?: string;
}

interface JobResult {
  job_id: string;
  status: "queued" | "running" | "done" | "error";
  steps: StepResult[];
  final?: { reactions: number; metabolites: number; genes: number };
  initial?: { reactions: number; metabolites: number; genes: number };
  reduction_pct?: number;
  duration_s?: number;
  ci_factor?: number;
  error?: string;
}

// ── Static pipeline descriptions ──────────────────────────────────────────────

interface PipelineStep {
  step: number;
  name: string;
  icon: string;
  method: string;
  desc: string;
  color: string;
  paperRxns?: number;
  tasks?: number;
  note?: string;
}

const PIPELINE: PipelineStep[] = [
  {
    step: 0, name: "Model Setup", icon: "⚙", method: "Manual curation", color: "#475569",
    desc: "Curate iCHO1766 for CHO-DG44 auxotrophies (Arg, Cys, Pro, Lys). Add fumarate/fucose exchange reactions missing from the original model. Set exchange bounds from 57 measured metabolites across all 12 cultures.",
    paperRxns: 6663, tasks: 155,
  },
  {
    step: 1, name: "Infeasibility Resolution", icon: "⚡", method: "Slack LP  (Eq. 3)", color: "#1d6fa5",
    desc: "For each time point, a slack-variable LP checks feasibility under imposed MetRaC bounds. If infeasible, minimal slack δ⁺/δ⁻ is added only to exchange reactions. Persistent infeasibilities get auxiliary demand reactions set to the minimum slack needed across all time points.",
    paperRxns: 4400, tasks: 136,
    note: "Network reduced by >1/3. No data-caused infeasibilities at 95% CI.",
  },
  {
    step: 2, name: "Exchange Pruning", icon: "🔍", method: "MILP  (Eq. 4)", color: "#2c9c56",
    desc: "Binary MILP (min Σwᵢ) identifies the smallest subset of exchange reactions needed for a feasible flux distribution at every time point. The MILP is solved separately at each time point and the union of active exchanges is retained.",
    tasks: 131,
    note: "Converges to 37 essential exchanges = 29 core metabolites + product + growth + 7 cofactors.",
  },
  {
    step: 3, name: "Transport Cleanup", icon: "🔄", method: "Biological rules", color: "#a05ca0",
    desc: "Each retained exchange typically has multiple possible transport mechanisms. Preference hierarchy: passive > proton-coupled > sodium-coupled. Redundant transport routes are removed.",
    note: "Eliminates MILP alternate-optima by fixing transport mechanism biologically.",
  },
  {
    step: 4, name: "pFBA Trimming", icon: "✂", method: "pFBA", color: "#e07b3c",
    desc: "Parsimonious FBA (minimise Σv² subject to max objective) is run at every experimental time point. Reactions that carry zero flux in every parsimonious solution are permanently removed.",
    paperRxns: 860, tasks: 105,
  },
  {
    step: 5, name: "Loop Removal", icon: "⟲", method: "Loopless FBA", color: "#c45252",
    desc: "Standard FBA vs. loopless FBA comparison across all time points and objectives (biomass + IgG1). Reactions active only under standard FBA are flagged and removed.",
    paperRxns: 860, tasks: 105,
    note: "Structure unchanged — prior steps had already eliminated loop artifacts.",
  },
];

const CI_DATA = [
  { ci: "68%",  reactions: 770,  demandRxns: 10, tag: "Needs demand rxns" },
  { ci: "95%",  reactions: 860,  demandRxns: 0,  tag: "Optimal ★" },
  { ci: "99%",  reactions: 1040, demandRxns: 0,  tag: "Wider bounds" },
  { ci: "100%", reactions: 1260, demandRxns: 0,  tag: "Exact mean" },
];

const SUBSYSTEM_DATA = [
  { name: "Exchange / Transport", pct: 92, color: "#1d6fa5" },
  { name: "Lipid Metabolism",     pct: 90, color: "#e07b3c" },
  { name: "Amino Acids",          pct: 73, color: "#2c9c56" },
  { name: "Nucleotides",          pct: 62, color: "#a05ca0" },
  { name: "Energy Metabolism",    pct: 34, color: "#c45252" },
];

const METHOD_CMP = [
  { method: "MetRaC 95% CI", reactions: 860, tasks: "105/155", demandRxns: "None ✓", best: true },
  { method: "Constant rate (α = 0.3)", reactions: 1075, tasks: "~100/155", demandRxns: "Several ✗", best: false },
  { method: "GIMME (transcriptomics)", reactions: 1200, tasks: "<105/155", demandRxns: "Many ✗", best: false },
];

const CORE_METS: { name: string; aux: boolean }[] = [
  "Acetate", "Alanine", "Asparagine", "Aspartate", "Citrate",
  "Cystine", "Formate", "Fumarate", "Glucose", "Glutamine",
  "Glutamate", "Glycine", "Histidine", "Isoleucine", "Lactate",
  "Leucine", "Methionine", "Ammonium", "Phenylalanine",
  "Serine", "Succinate", "Threonine", "Tryptophan", "Tyrosine", "Valine",
].map(n => ({ name: n, aux: false })).concat([
  { name: "Arginine", aux: true }, { name: "Cysteine", aux: true },
  { name: "Lysine", aux: true },   { name: "Proline",  aux: true },
]).sort((a, b) => a.name.localeCompare(b.name));

const COFACTORS = ["CO₂", "H⁺", "H₂O", "HCO₃⁻", "O₂", "SO₄²⁻", "Pᵢ"];

// ── Live Run Panel ─────────────────────────────────────────────────────────────

interface LivePanelProps {
  modelInfo: ModelInfo | null;
  modelLoading: boolean;
  onLoadModel: () => void;
}

function LivePanel({ modelInfo, modelLoading, onLoadModel }: LivePanelProps) {
  const [ciLevel, setCiLevel] = useState<number>(95);
  const [maxStep, setMaxStep] = useState<number>(5);
  const [job, setJob] = useState<JobResult | null>(null);
  const [running, setRunning] = useState(false);
  const [elapsed, setElapsed] = useState(0);
  const pollRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const timerRef = useRef<ReturnType<typeof setInterval> | null>(null);

  const ciToFactor = (ci: number) => ci >= 99 ? 0.0 : ci >= 95 ? 0.10 : ci >= 90 ? 0.18 : 0.30;

  const stopPolling = useCallback(() => {
    if (pollRef.current) { clearInterval(pollRef.current); pollRef.current = null; }
    if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null; }
  }, []);

  const pollJob = useCallback(async (jobId: string) => {
    try {
      const r = await fetch(`/gem/job/${jobId}`);
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const data: JobResult = await r.json();
      setJob(data);
      if (data.status === "done" || data.status === "error") {
        stopPolling();
        setRunning(false);
      }
    } catch (e) {
      console.error("Poll error", e);
    }
  }, [stopPolling]);

  const runPipeline = useCallback(async () => {
    stopPolling();
    setRunning(true);
    setJob(null);
    setElapsed(0);
    const t0 = Date.now();
    timerRef.current = setInterval(() => setElapsed(Math.round((Date.now() - t0) / 1000)), 500);

    try {
      const r = await fetch("/gem/run-pipeline", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ ci_factor: ciToFactor(ciLevel), max_step: maxStep, quick_mode: true }),
      });
      if (!r.ok) { const t = await r.text(); throw new Error(`HTTP ${r.status}: ${t}`); }
      const { job_id } = await r.json();
      setJob({ job_id, status: "queued", steps: [] });
      // Poll every 1.5s
      pollRef.current = setInterval(() => pollJob(job_id), 1500);
    } catch (e) {
      setJob({ job_id: "", status: "error", steps: [], error: String(e) });
      setRunning(false);
      stopPolling();
    }
  }, [ciLevel, maxStep, pollJob, stopPolling]);

  useEffect(() => () => stopPolling(), [stopPolling]);

  const compression = job?.reduction_pct;
  const isActive = running || (job?.status === "queued") || (job?.status === "running");

  return (
    <div className="gem-live-panel">
      <div className="gem-live-header">
        <span className="gem-live-badge">⚡ LIVE</span>
        <span className="gem-live-title">Run Pipeline on iCHOv1 — FastAPI + COBRApy Backend</span>
      </div>

      {/* Model info strip */}
      {modelLoading && <div className="gem-live-loading">Loading iCHOv1 model info…</div>}
      {modelInfo && (
        <div className="gem-model-strip">
          <div className="gem-model-chip">
            <span className="gem-model-val">{modelInfo.reactions.toLocaleString()}</span>
            <span className="gem-model-lab">reactions</span>
          </div>
          <div className="gem-model-chip">
            <span className="gem-model-val">{modelInfo.metabolites.toLocaleString()}</span>
            <span className="gem-model-lab">metabolites</span>
          </div>
          <div className="gem-model-chip">
            <span className="gem-model-val">{modelInfo.genes.toLocaleString()}</span>
            <span className="gem-model-lab">genes</span>
          </div>
          <div className="gem-model-chip">
            <span className="gem-model-val">{modelInfo.exchange_count}</span>
            <span className="gem-model-lab">exchanges</span>
          </div>
          <div className="gem-model-chip">
            <span className="gem-model-val">{modelInfo.core_exchanges_present.length}</span>
            <span className="gem-model-lab">core met. found</span>
          </div>
          <div className="gem-model-id">{modelInfo.model_id}</div>
        </div>
      )}
      {!modelInfo && !modelLoading && (
        <button className="gem-load-btn" onClick={onLoadModel}>Load model info</button>
      )}

      {/* Controls */}
      <div className="gem-controls-row">
        <div className="gem-control-group">
          <label className="gem-control-label">CI Level</label>
          <select className="gem-select" value={ciLevel} onChange={e => setCiLevel(Number(e.target.value))} disabled={isActive}>
            <option value={68}>68% CI</option>
            <option value={90}>90% CI</option>
            <option value={95}>95% CI (paper default)</option>
            <option value={99}>99% CI</option>
          </select>
        </div>
        <div className="gem-control-group">
          <label className="gem-control-label">Run through Step</label>
          <select className="gem-select" value={maxStep} onChange={e => setMaxStep(Number(e.target.value))} disabled={isActive}>
            {PIPELINE.map(s => (
              <option key={s.step} value={s.step}>Step {s.step} — {s.name}</option>
            ))}
          </select>
        </div>
        <button
          className={`gem-run-btn${isActive ? " gem-run-running" : ""}`}
          onClick={runPipeline}
          disabled={isActive}>
          {isActive ? `Running… ${elapsed}s` : job?.status === "done" ? "▶ Run Again" : "▶ Run Pipeline"}
        </button>
      </div>

      {isActive && (
        <div className="gem-progress-bar"><div className="gem-progress-fill" /></div>
      )}

      {/* Live step progress during run */}
      {job && job.steps.length > 0 && job.status !== "done" && (
        <div className="gem-live-steps">
          {job.steps.map(s => (
            <div key={s.step} className="gem-live-step-row">
              <span className="gem-step-pill" style={{ background: PIPELINE[s.step]?.color }}>S{s.step}</span>
              <span className="gem-live-step-name">{s.name}</span>
              <span className="gem-live-step-rxns">{s.reactions.toLocaleString()} rxns</span>
              <span className="gem-live-step-time">{s.duration_s?.toFixed(2)}s</span>
            </div>
          ))}
          {job.status === "running" && (
            <div className="gem-live-step-row gem-live-step-pending">
              <span className="gem-step-pill" style={{ background: "#94a3b8" }}>…</span>
              <span className="gem-live-step-name">Running next step…</span>
            </div>
          )}
        </div>
      )}

      {job?.status === "error" && (
        <div className="gem-live-error"><strong>Error:</strong> {job.error}</div>
      )}

      {/* Final results */}
      {job?.status === "done" && job.final && (
        <div className="gem-results">
          <div className="gem-results-kpis">
            <div className="gem-res-kpi">
              <span className="gem-res-kpi-val gem-kpi-green">{job.final.reactions.toLocaleString()}</span>
              <span className="gem-res-kpi-lab">reactions</span>
            </div>
            <div className="gem-res-kpi">
              <span className="gem-res-kpi-val">{job.final.metabolites.toLocaleString()}</span>
              <span className="gem-res-kpi-lab">metabolites</span>
            </div>
            <div className="gem-res-kpi">
              <span className="gem-res-kpi-val">{job.final.genes.toLocaleString()}</span>
              <span className="gem-res-kpi-lab">genes</span>
            </div>
            {compression !== undefined && (
              <div className="gem-res-kpi">
                <span className="gem-res-kpi-val gem-kpi-red">−{compression}%</span>
                <span className="gem-res-kpi-lab">compression</span>
              </div>
            )}
            <div className="gem-res-kpi">
              <span className="gem-res-kpi-val">{job.duration_s?.toFixed(1)}s</span>
              <span className="gem-res-kpi-lab">wall time</span>
            </div>
          </div>

          {/* Step-by-step table */}
          <div className="gem-step-table-wrap">
            <table className="gem-step-table">
              <thead>
                <tr>
                  <th>Step</th><th>Name</th>
                  <th className="gem-td-num">Reactions</th>
                  <th className="gem-td-num">Metabolites</th>
                  <th className="gem-td-num">Genes</th>
                  <th className="gem-td-num">Time (s)</th>
                  <th>Notes</th>
                </tr>
              </thead>
              <tbody>
                {job.steps.map(s => {
                  const paper = PIPELINE[s.step]?.paperRxns;
                  const diff = paper !== undefined ? s.reactions - paper : null;
                  return (
                    <tr key={s.step}>
                      <td><span className="gem-step-pill" style={{ background: PIPELINE[s.step]?.color }}>{s.step}</span></td>
                      <td>{s.name}</td>
                      <td className="gem-td-num">
                        {s.reactions.toLocaleString()}
                        {diff !== null && (
                          <span className={`gem-diff ${diff > 0 ? "gem-diff-pos" : diff < 0 ? "gem-diff-neg" : "gem-diff-zero"}`}>
                            {diff > 0 ? `+${diff}` : diff < 0 ? `${diff}` : "="}
                          </span>
                        )}
                      </td>
                      <td className="gem-td-num">{s.metabolites.toLocaleString()}</td>
                      <td className="gem-td-num">{s.genes.toLocaleString()}</td>
                      <td className="gem-td-num">{s.duration_s !== undefined ? s.duration_s.toFixed(2) : "—"}</td>
                      <td className="gem-step-note-cell">
                        {s.demand_reactions !== undefined && s.demand_reactions > 0 && (
                          <span className="gem-note-chip gem-note-warn">{s.demand_reactions} demand rxns</span>
                        )}
                        {s.infeasible_timepoints !== undefined && s.infeasible_timepoints > 0 && (
                          <span className="gem-note-chip gem-note-warn">{s.infeasible_timepoints} infeasible tps</span>
                        )}
                        {s.exchanges_kept !== undefined && (
                          <span className="gem-note-chip">{s.exchanges_kept} exchanges</span>
                        )}
                        {s.note && <span className="gem-note-chip gem-note-info">{s.note}</span>}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          {/* Live vs paper waterfall */}
          <div className="gem-waterfall-wrap">
            <div className="gem-waterfall-label">Reactions per step — live result vs. paper target</div>
            <ResponsiveContainer width="100%" height={160}>
              <BarChart
                data={job.steps.map(s => ({
                  name: `S${s.step}`,
                  live: s.reactions,
                  paper: PIPELINE[s.step]?.paperRxns,
                }))}
                margin={{ top: 4, right: 8, left: 0, bottom: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(0,0,0,0.07)" />
                <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                <YAxis tick={{ fontSize: 10 }} width={52} />
                <Tooltip contentStyle={{ fontSize: 10 }}
                  formatter={(v: number, name: string) => [v?.toLocaleString(), name === "live" ? "Live (this run)" : "Paper"]} />
                <Bar dataKey="paper" name="paper" fill="#cbd5e1" radius={[2, 2, 0, 0]} />
                <Bar dataKey="live"  name="live"  fill="#1d4ed8" radius={[2, 2, 0, 0]} />
              </BarChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </div>
  );
}

// ── Page ───────────────────────────────────────────────────────────────────────

export default function GEMReductionPage() {
  const [activeStep, setActiveStep] = useState<number | null>(null);
  const [modelInfo, setModelInfo] = useState<ModelInfo | null>(null);
  const [modelLoading, setModelLoading] = useState(false);

  const loadModelInfo = useCallback(async () => {
    setModelLoading(true);
    try {
      const r = await fetch("/gem/model-info");
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      setModelInfo(await r.json());
    } catch (e) {
      console.error("Failed to load model info", e);
    } finally {
      setModelLoading(false);
    }
  }, []);

  useEffect(() => { loadModelInfo(); }, [loadModelInfo]);

  const sel = activeStep !== null ? PIPELINE[activeStep] : null;

  return (
    <div className="gem-page">

      {/* Header */}
      <div className="gem-header">
        <div>
          <h1 className="gem-title">GEM Reduction Pipeline</h1>
          <p className="gem-subtitle">
            Uncertainty-aware iCHO1766 → compact model for digital twin integration ·{" "}
            <a href="https://doi.org/10.1038/s41540-026-00704-4" target="_blank" rel="noopener" className="gem-paper-link">
              Antonakoudis & Richelle (2026) npj Syst. Biol. Appl.
            </a>
          </p>
        </div>
        <div className="gem-kpi-row">
          <div className="gem-kpi"><span className="gem-kpi-val">6,663</span><span className="gem-kpi-lab">original reactions</span></div>
          <span className="gem-kpi-arrow">→</span>
          <div className="gem-kpi"><span className="gem-kpi-val gem-kpi-green">860</span><span className="gem-kpi-lab">reduced reactions</span></div>
          <div className="gem-kpi"><span className="gem-kpi-val gem-kpi-red">−87%</span><span className="gem-kpi-lab">compression</span></div>
          <div className="gem-kpi">
            <span className="gem-kpi-val">105<span style={{ fontWeight: 400, fontSize: "0.75rem" }}>/155</span></span>
            <span className="gem-kpi-lab">tasks retained</span>
          </div>
        </div>
      </div>

      {/* Algorithm box */}
      <div className="gem-algo-box">
        <strong>Pipeline:</strong>&nbsp;
        (0) Curate iCHO1766 → (1) Slack LP resolves infeasibilities → (2) MILP selects 37 essential exchanges →
        (3) Transport deduplication → (4) pFBA zeros out unused reactions → (5) Loopless FBA removes cycles.
        MetRaC 95% CI bounds constrain all steps.
        &nbsp;<span className="gem-algo-note">
          Implemented in COBRApy with GLPK LP/MILP solver. Applied to 12 CHO DG44 fed-batch cultures.
        </span>
      </div>

      {/* LIVE run panel */}
      <LivePanel modelInfo={modelInfo} modelLoading={modelLoading} onLoadModel={loadModelInfo} />

      {/* Pipeline steps */}
      <div className="gem-section-label">5-Step Pipeline — click any step for details</div>
      <div className="gem-pipeline-track">
        {PIPELINE.map((step, i) => (
          <React.Fragment key={step.step}>
            <div
              className={`gem-step-card${activeStep === i ? " gem-step-active" : ""}`}
              style={{ "--sc": step.color } as React.CSSProperties}
              onClick={() => setActiveStep(activeStep === i ? null : i)}>
              <div className="gem-step-icon">{step.icon}</div>
              <div className="gem-step-num">Step {step.step}</div>
              <div className="gem-step-name">{step.name}</div>
              <div className="gem-step-method">{step.method}</div>
              {step.paperRxns !== undefined && <div className="gem-step-rxn">{step.paperRxns.toLocaleString()} rxns</div>}
              {step.tasks !== undefined && <div className="gem-step-tasks">{step.tasks}/155 tasks</div>}
            </div>
            {i < PIPELINE.length - 1 && <div className="gem-step-arrow">▶</div>}
          </React.Fragment>
        ))}
      </div>

      {sel && (
        <div className="gem-step-detail" style={{ borderLeftColor: sel.color }}>
          <div className="gem-step-detail-hd">
            <strong>Step {sel.step} — {sel.name}</strong>
            <span className="gem-method-badge" style={{ background: sel.color }}>{sel.method}</span>
          </div>
          <p className="gem-step-detail-body">{sel.desc}</p>
          {sel.note && <p className="gem-step-detail-note">💡 {sel.note}</p>}
          {(sel.paperRxns !== undefined || sel.tasks !== undefined) && (
            <div className="gem-step-chips">
              {sel.paperRxns !== undefined && <span className="gem-chip">{sel.paperRxns.toLocaleString()} reactions (paper)</span>}
              {sel.tasks !== undefined && <span className="gem-chip">{sel.tasks}/155 tasks</span>}
            </div>
          )}
        </div>
      )}

      {/* Chart row */}
      <div className="gem-charts-row">
        <div className="gem-chart-card">
          <h3 className="gem-chart-title">CI Threshold vs. Network Size</h3>
          <p className="gem-chart-sub">MetRaC bounds at four CI levels. 95% CI gives the smallest model that is fully feasible without extra demand reactions.</p>
          <ResponsiveContainer width="100%" height={180}>
            <BarChart data={CI_DATA} margin={{ top: 4, right: 8, left: 0, bottom: 0 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(0,0,0,0.07)" />
              <XAxis dataKey="ci" tick={{ fontSize: 11 }} />
              <YAxis tick={{ fontSize: 10 }} width={52}
                label={{ value: "Reactions", angle: -90, position: "insideLeft", fontSize: 9, dx: 10 }} />
              <Tooltip
                formatter={(v: number, _n: string, p: { payload: typeof CI_DATA[0] }) => [
                  `${v.toLocaleString()} (${p.payload.demandRxns > 0 ? p.payload.demandRxns + " demand rxns" : "no demand rxns"})`,
                  "Reactions",
                ]}
                contentStyle={{ fontSize: 10 }} />
              <Bar dataKey="reactions" radius={[3, 3, 0, 0]}>
                {CI_DATA.map((d, i) => (
                  <Cell key={i} fill={d.ci === "95%" ? "#1d4ed8" : d.demandRxns > 0 ? "#f59e0b" : "#94a3b8"} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
          <div className="gem-ci-legend">
            {CI_DATA.map(d => (
              <div key={d.ci} className={`gem-ci-chip${d.ci === "95%" ? " gem-ci-opt" : d.demandRxns > 0 ? " gem-ci-warn" : ""}`}>
                <strong>{d.ci}</strong> {d.tag}
              </div>
            ))}
          </div>
        </div>

        <div className="gem-chart-card">
          <h3 className="gem-chart-title">Subsystem Compression</h3>
          <p className="gem-chart-sub">% of reactions removed per metabolic subsystem at 95% CI.</p>
          <div className="gem-subsys-list">
            {SUBSYSTEM_DATA.map(d => (
              <div key={d.name} className="gem-subsys-row">
                <span className="gem-subsys-name">{d.name}</span>
                <div className="gem-subsys-track"><div className="gem-subsys-fill" style={{ width: `${d.pct}%`, background: d.color }} /></div>
                <span className="gem-subsys-pct">−{d.pct}%</span>
              </div>
            ))}
          </div>
        </div>

        <div className="gem-chart-card">
          <h3 className="gem-chart-title">Method Comparison</h3>
          <p className="gem-chart-sub">MetRaC uncertainty-aware reduction vs. constant-rate and GIMME baselines.</p>
          <table className="gem-cmp-table">
            <thead><tr><th>Method</th><th>Rxns</th><th>Tasks</th><th>Demand</th></tr></thead>
            <tbody>
              {METHOD_CMP.map(m => (
                <tr key={m.method} className={m.best ? "gem-cmp-best" : ""}>
                  <td>{m.method}{m.best && <span className="gem-best-badge">★ Best</span>}</td>
                  <td className="gem-td-num">{m.reactions.toLocaleString()}</td>
                  <td className="gem-td-num">{m.tasks}</td>
                  <td className="gem-td-num">{m.demandRxns}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Lower row */}
      <div className="gem-lower-row">
        <div className="gem-mets-card">
          <h3 className="gem-chart-title">37 Essential Exchanges — Step 2 MILP Result</h3>
          <p className="gem-chart-sub">★ = CHO-DG44 auxotrophic amino acid (pathway blocked).</p>
          <div className="gem-mets-block">
            <div className="gem-mets-label">29 Core Metabolites</div>
            <div className="gem-mets-grid">
              {CORE_METS.map(m => (
                <div key={m.name} className={`gem-met-chip${m.aux ? " gem-met-aux" : ""}`}>
                  {m.name}{m.aux && <span className="gem-aux-star">★</span>}
                </div>
              ))}
            </div>
          </div>
          <div className="gem-mets-block" style={{ marginTop: "0.6rem" }}>
            <div className="gem-mets-label">Essential Cofactors + Product/Growth</div>
            <div className="gem-mets-grid">
              {COFACTORS.map(m => <div key={m} className="gem-met-chip gem-met-cofactor">{m}</div>)}
              <div className="gem-met-chip gem-met-product">IgG1 mAb</div>
              <div className="gem-met-chip gem-met-product">μ_eff</div>
            </div>
          </div>
        </div>

        <div className="gem-eq-card">
          <h3 className="gem-chart-title">Key Optimisation Formulations</h3>
          <div className="gem-eq-block">
            <div className="gem-eq-label">Step 1 — Infeasibility slack LP (Eq. 3)</div>
            <pre className="gem-eq-pre">{`min   Σᵢ (δᵢ⁻ + δᵢ⁺)
s.t.  S · v  = 0
      vᵢ ≥ LBᵢ − δᵢ⁻    ∀i
      vᵢ ≤ UBᵢ + δᵢ⁺    ∀i
      δᵢ⁻, δᵢ⁺ ≥ 0
      δ  = 0  if i ∉ exchange reactions`}</pre>
          </div>
          <div className="gem-eq-block">
            <div className="gem-eq-label">Step 2 — Exchange pruning MILP (Eq. 4)</div>
            <pre className="gem-eq-pre">{`min   Σᵢ wᵢ               (wᵢ ∈ {0,1})
s.t.  S · v = 0
      vᵢ ≤ wᵢ · UBᵢ       ∀i ∈ removable
      vᵢ ≥ wᵢ · LBᵢ + ε   ∀i ∈ removable`}</pre>
          </div>
          <div className="gem-eq-note">
            <strong>MetRaC vs. constant rate:</strong> MetRaC models each rate as a linear combination
            of logistic basis functions with posterior estimated via nested sampling — giving time-resolved
            CI bounds rather than point estimates with arbitrary α bounds.
          </div>
        </div>
      </div>

      {/* PC-dFBA link */}
      <div className="gem-link-box">
        <span className="gem-link-icon">⬡</span>
        <div>
          <strong>Connection to PC-dFBA:</strong> The reduced iCHO1766 model (860 reactions) produced
          here is the intended input to the PC-dFBA framework (Richelle et al. 2025). Flux distributions
          from this model are projected onto principal metabolic coordinates via PCA, then used to train
          the PC-dFBA neural network.
        </div>
      </div>

    </div>
  );
}
