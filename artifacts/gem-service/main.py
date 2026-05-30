"""
GEM Reduction Service — FastAPI backend
Antonakoudis & Richelle (2026) 5-step pipeline on iCHOv1 (Hefzi et al. 2016)

Job-queue pattern: POST /gem/run-pipeline returns a job_id immediately.
Client polls GET /gem/job/{job_id} until status == "done" or "error".
This avoids proxy / browser timeout issues for the ~25s pipeline run.
"""

from __future__ import annotations

import asyncio
import logging
import pickle
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock
from typing import Any

import numpy as np
import cobra
from cobra import Model, Reaction
from cobra.flux_analysis import pfba
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("gem-service")

# ── App ────────────────────────────────────────────────────────────────────────
app = FastAPI(title="GEM Reduction Service", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Thread pool for CPU-bound COBRApy work
_executor = ThreadPoolExecutor(max_workers=2)

# ── In-memory job store ────────────────────────────────────────────────────────
_jobs: dict[str, dict] = {}
_jobs_lock = Lock()

# ── Model loading ──────────────────────────────────────────────────────────────
MODEL_PATH = Path(__file__).parent / "models" / "iCHOv1.json"
_base_model: Model | None = None
_base_model_pickle: bytes | None = None

# Cached reduced model (stored after a successful GEM reduction pipeline run)
_reduced_model_pickle: bytes | None = None
_reduced_model_lock = Lock()


def _store_reduced_model(m: Model) -> None:
    global _reduced_model_pickle
    with _reduced_model_lock:
        _reduced_model_pickle = pickle.dumps(m)
        log.info("Cached reduced model (%d rxns)", len(m.reactions))


def _get_reduced_model() -> Model | None:
    with _reduced_model_lock:
        if _reduced_model_pickle is None:
            return None
        return pickle.loads(_reduced_model_pickle)


def get_base_model() -> Model:
    global _base_model, _base_model_pickle
    if _base_model is None:
        log.info("Loading iCHOv1 from %s …", MODEL_PATH)
        t0 = time.time()
        _base_model = cobra.io.load_json_model(str(MODEL_PATH))
        log.info("Loaded iCHOv1: %d rxns, %d mets in %.1fs",
                 len(_base_model.reactions), len(_base_model.metabolites), time.time() - t0)
        t1 = time.time()
        _base_model_pickle = pickle.dumps(_base_model)
        log.info("Pickled iCHOv1 in %.2fs (%d bytes)", time.time() - t1, len(_base_model_pickle))
    return pickle.loads(_base_model_pickle)


# ── CHO-DG44 auxotrophies ──────────────────────────────────────────────────────
CHO_DG44_BLOCKED = [
    "ARGSS", "CYSTGL", "CYSTS", "GLU5Km", "ORNTArm", "BTND1",
]

CORE_EXCHANGE_PREFIXES = [
    "EX_ac_e", "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e",
    "EX_cit_e", "EX_cys__L_e", "EX_cystin_e", "EX_for_e", "EX_fum_e",
    "EX_glc__D_e", "EX_gln__L_e", "EX_glu__L_e", "EX_gly_e", "EX_his__L_e",
    "EX_ile__L_e", "EX_lac__L_e", "EX_leu__L_e", "EX_lys__L_e", "EX_met__L_e",
    "EX_nh4_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e", "EX_succ_e",
    "EX_thr__L_e", "EX_trp__L_e", "EX_tyr__L_e", "EX_val__L_e",
    "EX_h_e", "EX_h2o_e", "EX_pi_e", "EX_co2_e", "EX_o2_e", "EX_so4_e", "EX_hco3_e",
]


def model_summary(m: Model) -> dict:
    return {
        "reactions": len(m.reactions),
        "metabolites": len(m.metabolites),
        "genes": len(m.genes),
        "exchange_reactions": len(m.exchanges),
    }


# ── Pipeline steps ─────────────────────────────────────────────────────────────

def step0_setup(m: Model) -> tuple[Model, dict]:
    t0 = time.time()
    blocked = []
    for rxn_id in CHO_DG44_BLOCKED:
        if rxn_id in m.reactions:
            m.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
            blocked.append(rxn_id)
    core_present = [r.id for r in m.exchanges
                    if any(r.id.startswith(p) for p in CORE_EXCHANGE_PREFIXES)]
    stats = model_summary(m)
    return m, {
        "step": 0, "name": "Model Setup",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "blocked_reactions": blocked, "exchanges_kept": len(m.exchanges),
        "core_exchanges_found": len(core_present),
        "duration_s": round(time.time() - t0, 2),
        "note": f"{len(blocked)} auxotrophy reactions blocked",
    }


def step1_infeasibility(m: Model, ci_factor: float = 0.1) -> tuple[Model, dict]:
    t0 = time.time()
    demand_added = []
    sol = m.optimize()
    if sol.status == "infeasible":
        for rxn in m.exchanges:
            demand = Reaction(f"DEMAND_{rxn.id}")
            demand.lower_bound = 0
            demand.upper_bound = 1000
            met = list(rxn.metabolites.keys())[0]
            demand.add_metabolites({met: -1})
            m.add_reactions([demand])
            demand_added.append(rxn.id)

    infeasible_count = 0
    for rxn in list(m.exchanges)[:50]:
        orig_lb, orig_ub = rxn.bounds
        new_lb = orig_lb * (1 - ci_factor) if orig_lb < 0 else orig_lb * (1 + ci_factor)
        new_ub = orig_ub * (1 - ci_factor) if orig_ub > 0 else orig_ub * (1 + ci_factor)
        with m:
            rxn.bounds = (new_lb, new_ub)
            s = m.optimize()
            if s.status == "infeasible":
                infeasible_count += 1

    dead_rxns = [r for r in m.reactions if r.bounds == (0.0, 0.0)
                 and r.id not in {r2.id for r2 in m.exchanges}]
    m.remove_reactions(dead_rxns[:200], remove_orphans=True)

    stats = model_summary(m)
    return m, {
        "step": 1, "name": "Infeasibility Resolution",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "infeasible_timepoints": infeasible_count, "demand_reactions": len(demand_added),
        "exchanges_kept": len(m.exchanges), "duration_s": round(time.time() - t0, 2),
        "note": f"CI factor {ci_factor:.0%}; {len(dead_rxns[:200])} zero-flux rxns removed",
    }


def step2_exchange_milp(m: Model, quick: bool = True) -> tuple[Model, dict]:
    t0 = time.time()
    core_ids = set(CORE_EXCHANGE_PREFIXES)
    candidate_removals, always_keep = [], []
    for rxn in m.exchanges:
        if any(rxn.id.startswith(p) for p in core_ids):
            always_keep.append(rxn.id)
        else:
            candidate_removals.append(rxn)

    removed, kept_necessary = [], []
    max_test = min(len(candidate_removals), 100) if quick else len(candidate_removals)
    for rxn in candidate_removals[:max_test]:
        with m:
            rxn.bounds = (0.0, 0.0)
            if m.optimize().status == "optimal":
                removed.append(rxn.id)
            else:
                kept_necessary.append(rxn.id)
    for rxn in candidate_removals[max_test:]:
        kept_necessary.append(rxn.id)

    m.remove_reactions([m.reactions.get_by_id(rid) for rid in removed if rid in m.reactions],
                       remove_orphans=True)
    stats = model_summary(m)
    return m, {
        "step": 2, "name": "Exchange Pruning (MILP)",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "exchanges_kept": len(always_keep) + len(kept_necessary), "exchanges_removed": len(removed),
        "duration_s": round(time.time() - t0, 2),
        "note": f"{'Quick: ' if quick else ''}{len(removed)} exchanges removed, {len(always_keep)+len(kept_necessary)} kept",
    }


def step3_transport(m: Model) -> tuple[Model, dict]:
    t0 = time.time()
    removed = []
    met_to_transport: dict[str, list[Reaction]] = {}
    for rxn in list(m.reactions):
        compartments = {met.compartment for met in rxn.metabolites}
        if len(compartments) > 1 and "e" in compartments:
            for met in rxn.metabolites:
                if met.compartment == "e":
                    met_to_transport.setdefault(met.id, []).append(rxn)

    for met_id, rxns in met_to_transport.items():
        if len(rxns) <= 1:
            continue
        def score(rxn: Reaction) -> int:
            ids = {met.id for met in rxn.metabolites}
            has_h = any("h_c" in mid or mid == "h_e" for mid in ids)
            has_na = any("na1" in mid.lower() for mid in ids)
            return 0 if not has_h and not has_na else (1 if has_h and not has_na else 2)
        for rxn in sorted(rxns, key=score)[1:]:
            if rxn.id in m.reactions:
                removed.append(rxn.id)

    m.remove_reactions([m.reactions.get_by_id(rid) for rid in removed if rid in m.reactions],
                       remove_orphans=True)
    stats = model_summary(m)
    return m, {
        "step": 3, "name": "Transport Cleanup",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "exchanges_kept": len(m.exchanges), "duration_s": round(time.time() - t0, 2),
        "note": f"{len(removed)} redundant transport rxns removed",
    }


def step4_pfba(m: Model) -> tuple[Model, dict]:
    t0 = time.time()
    try:
        sol = pfba(m)
        status = sol.status
    except Exception as exc:
        log.warning("pFBA failed: %s — using standard FBA", exc)
        sol = m.optimize()
        status = sol.status + "_fallback"

    to_remove_ids = []
    if sol.fluxes is not None:
        protected = {r.id for r in m.exchanges} | {"BIOMASS_cho", "BIOMASS_CHO", "Growth"}
        zero_flux = [rxn.id for rxn in m.reactions if abs(sol.fluxes.get(rxn.id, 0.0)) < 1e-9]
        to_remove_ids = [rid for rid in zero_flux if rid not in protected]
        m.remove_reactions([m.reactions.get_by_id(rid) for rid in to_remove_ids if rid in m.reactions],
                           remove_orphans=True)
    stats = model_summary(m)
    return m, {
        "step": 4, "name": "pFBA Trimming",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "exchanges_kept": len(m.exchanges), "duration_s": round(time.time() - t0, 2),
        "note": f"pFBA {status}; {len(to_remove_ids)} zero-flux rxns removed",
    }


def step5_loops(m: Model) -> tuple[Model, dict]:
    t0 = time.time()
    removed_ids = []
    try:
        from cobra.flux_analysis import loopless_solution
        std_sol = m.optimize()
        loop_sol = loopless_solution(m)
        if std_sol.fluxes is not None and loop_sol.fluxes is not None:
            protected = {r.id for r in m.exchanges}
            for rxn in list(m.reactions):
                if (abs(std_sol.fluxes.get(rxn.id, 0.0)) > 1e-9
                        and abs(loop_sol.fluxes.get(rxn.id, 0.0)) < 1e-9
                        and rxn.id not in protected):
                    removed_ids.append(rxn.id)
            m.remove_reactions([m.reactions.get_by_id(rid) for rid in removed_ids if rid in m.reactions],
                               remove_orphans=True)
        loop_status = "ok"
    except Exception as exc:
        log.warning("Loopless FBA failed: %s", exc)
        loop_status = f"skipped ({exc})"

    stats = model_summary(m)
    return m, {
        "step": 5, "name": "Loop Removal",
        "reactions": stats["reactions"], "metabolites": stats["metabolites"], "genes": stats["genes"],
        "exchanges_kept": len(m.exchanges), "duration_s": round(time.time() - t0, 2),
        "note": f"Loopless FBA {loop_status}; {len(removed_ids)} cycle rxns removed",
    }


# ── Sync pipeline runner (runs in thread pool) ────────────────────────────────

def _run_pipeline_sync(job_id: str, ci_factor: float, max_step: int, quick_mode: bool):
    """Runs the full pipeline synchronously. Called from thread pool. Updates job store."""
    def _update(patch: dict):
        with _jobs_lock:
            _jobs[job_id].update(patch)

    t_total = time.time()
    _update({"status": "running", "steps": [], "started_at": t_total})

    try:
        m = get_base_model()
        initial_stats = model_summary(m)

        steps_results = []

        def do_step(fn, *args):
            result_m, r = fn(*args)
            steps_results.append(r)
            _update({"steps": list(steps_results)})
            return result_m

        m = do_step(step0_setup, m)
        if max_step >= 1:
            m = do_step(step1_infeasibility, m, ci_factor)
        if max_step >= 2:
            m = do_step(step2_exchange_milp, m, quick_mode)
        if max_step >= 3:
            m = do_step(step3_transport, m)
        if max_step >= 4:
            m = do_step(step4_pfba, m)
        if max_step >= 5:
            m = do_step(step5_loops, m)

        final = model_summary(m)
        reduction_pct = round(100 * (1 - final["reactions"] / initial_stats["reactions"]), 1)
        _store_reduced_model(m)  # cache for pcdfba sampling
        _update({
            "status": "done",
            "final": final,
            "initial": initial_stats,
            "reduction_pct": reduction_pct,
            "duration_s": round(time.time() - t_total, 1),
            "ci_factor": ci_factor,
        })

    except Exception as exc:
        log.exception("Pipeline job %s failed", job_id)
        _update({"status": "error", "error": str(exc), "duration_s": round(time.time() - t_total, 1)})


# ── PC-dFBA pFBA sampling ─────────────────────────────────────────────────────
# Sample pFBA optima across a physiological grid of Glc/Gln exchange rates,
# then compute PCA to get GEM-derived PC loadings for the PC-dFBA tab.

PCDFBA_GLC_GRID = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0]   # mmol/(gDW·h), Glc uptake
PCDFBA_GLN_GRID = [0.02, 0.05, 0.1, 0.2, 0.4]        # mmol/(gDW·h), Gln uptake


def _run_pcdfba_sampling_sync(job_id: str) -> None:
    """Sample pFBA over a Glc×Gln grid on the reduced model, run PCA, store result."""
    def _update(patch: dict):
        with _jobs_lock:
            _jobs[job_id].update(patch)

    _update({"status": "running", "started_at": time.time()})
    t0 = time.time()

    try:
        m = _get_reduced_model()
        model_source = "reduced"
        if m is None:
            log.info("No reduced model cached — using base model for pcdfba sampling")
            m = get_base_model()
            m, _ = step0_setup(m)
            model_source = "base (unreduced)"

        # Locate key exchange reactions
        glc_rxn = next(
            (r for r in m.exchanges if r.id.startswith("EX_glc")), None
        )
        gln_rxn = next(
            (r for r in m.exchanges if r.id.startswith("EX_gln")), None
        )
        if glc_rxn is None or gln_rxn is None:
            _update({"status": "error",
                     "error": "Could not find EX_glc / EX_gln exchange reactions"})
            return

        # Fixed reaction order — all non-exchange reactions in the model
        rxn_order = [r.id for r in m.reactions if r not in m.exchanges]

        flux_matrix: list[list[float]] = []
        conditions: list[dict] = []

        for q_glc in PCDFBA_GLC_GRID:
            for q_gln in PCDFBA_GLN_GRID:
                with m:
                    glc_rxn.lower_bound = -abs(q_glc)
                    gln_rxn.lower_bound = -abs(q_gln)
                    try:
                        sol = pfba(m)
                        if sol.status != "optimal":
                            continue
                    except Exception:
                        continue
                fluxes = [float(sol.fluxes.get(rid, 0.0)) for rid in rxn_order]
                flux_matrix.append(fluxes)
                conditions.append({"q_glc": q_glc, "q_gln": q_gln})

        if len(flux_matrix) < 4:
            _update({"status": "error",
                     "error": f"Only {len(flux_matrix)} feasible conditions — PCA requires ≥4"})
            return

        # PCA via SVD
        X = np.array(flux_matrix)          # (n_cond, n_rxns)
        X_mean = X.mean(axis=0)
        X_c = X - X_mean

        # Filter zero-variance reactions before SVD
        keep_mask = X_c.var(axis=0) > 1e-9
        X_f = X_c[:, keep_mask]
        rxn_order_f = [rid for rid, k in zip(rxn_order, keep_mask) if k]

        U, S, Vt = np.linalg.svd(X_f, full_matrices=False)
        scores = (U * S).tolist()             # (n_cond, min(n_cond, n_rxns_f))
        var_explained = (
            (S[:3] ** 2 / (S ** 2).sum()).tolist() if len(S) > 0 else []
        )

        def top_loadings(pc_idx: int, n_top: int = 20) -> list[dict]:
            if pc_idx >= len(Vt):
                return []
            loads = Vt[pc_idx]
            idx = np.argsort(np.abs(loads))[::-1][:n_top]
            return [
                {"rxn": rxn_order_f[int(i)], "loading": float(loads[i])}
                for i in idx
            ]

        _update({
            "status": "done",
            "n_conditions": len(flux_matrix),
            "n_reactions": len(rxn_order_f),
            "model_source": model_source,
            "var_explained": var_explained[:3],
            "scores": [s[:3] for s in scores],
            "conditions": conditions,
            "pc1_loadings": top_loadings(0),
            "pc2_loadings": top_loadings(1),
            "duration_s": round(time.time() - t0, 1),
        })

    except Exception as exc:
        log.exception("PC-dFBA sampling job %s failed", job_id)
        _update({"status": "error", "error": str(exc),
                 "duration_s": round(time.time() - t0, 1)})


# ── API routes ─────────────────────────────────────────────────────────────────

class RunPipelineRequest(BaseModel):
    ci_factor: float = 0.10
    max_step: int = 5
    quick_mode: bool = True


@app.get("/gem/healthz")
async def healthz():
    return {"status": "ok", "model_loaded": _base_model is not None}


@app.get("/gem/model-info")
async def model_info():
    try:
        loop = asyncio.get_event_loop()
        m = await loop.run_in_executor(_executor, get_base_model)
        return {
            "model_id": m.id,
            "reactions": len(m.reactions),
            "metabolites": len(m.metabolites),
            "genes": len(m.genes),
            "exchange_count": len(m.exchanges),
            "sample_exchanges": [r.id for r in m.exchanges[:50]],
            "core_exchanges_present": [
                r.id for r in m.exchanges
                if any(r.id.startswith(p) for p in CORE_EXCHANGE_PREFIXES)
            ],
        }
    except Exception as exc:
        log.exception("model-info failed")
        raise HTTPException(500, str(exc))


@app.post("/gem/run-pipeline")
async def start_pipeline(req: RunPipelineRequest):
    """Start the pipeline asynchronously. Returns job_id to poll."""
    job_id = str(uuid.uuid4())[:8]
    with _jobs_lock:
        _jobs[job_id] = {
            "job_id": job_id,
            "status": "queued",
            "steps": [],
            "ci_factor": req.ci_factor,
            "max_step": req.max_step,
        }
    # Fire off in thread pool — does not block the response
    loop = asyncio.get_event_loop()
    loop.run_in_executor(
        _executor,
        _run_pipeline_sync,
        job_id, req.ci_factor, req.max_step, req.quick_mode,
    )
    return {"job_id": job_id, "status": "queued"}


@app.get("/gem/job/{job_id}")
async def get_job(job_id: str):
    """Poll this endpoint for job status and incremental step results."""
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    return job


@app.delete("/gem/job/{job_id}")
async def delete_job(job_id: str):
    with _jobs_lock:
        _jobs.pop(job_id, None)
    return {"ok": True}


@app.post("/gem/pcdfba/run-sampling")
async def start_pcdfba_sampling():
    """Sample pFBA over Glc×Gln grid and compute PCA. Returns job_id to poll."""
    job_id = "pcd_" + str(uuid.uuid4())[:8]
    with _jobs_lock:
        _jobs[job_id] = {"job_id": job_id, "status": "queued"}
    loop = asyncio.get_event_loop()
    loop.run_in_executor(_executor, _run_pcdfba_sampling_sync, job_id)
    return {"job_id": job_id, "status": "queued"}


# ── Startup: pre-load model ────────────────────────────────────────────────────
@app.on_event("startup")
async def startup_event():
    log.info("Pre-loading iCHOv1 model …")
    try:
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(_executor, get_base_model)
        log.info("iCHOv1 pre-load complete")
    except Exception as exc:
        log.error("Failed to pre-load model: %s", exc)
