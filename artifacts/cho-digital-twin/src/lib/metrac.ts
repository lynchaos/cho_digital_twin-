
/**
 * MetRaC — Metabolic Rate Calculation  (§2.2 Richelle et al. 2025)
 *
 * Two estimation methods:
 *
 *   "kernel" (default, fast):
 *     1. Finite-difference dC/dt on noisy concentration measurements
 *     2. Bolus correction: intervals spanning a feed event are de-weighted (10× σ)
 *     3. q_i = -(dC_i/dt) / Xv
 *     4. Gaussian error propagation → σ_q
 *     5. Nadaraya-Watson kernel smoother → posterior mean + 95% CI
 *
 *   "gp" (Bayesian, proper):
 *     1. Fit SE-kernel Gaussian Process to raw concentration measurements
 *     2. Length-scale l optimised by marginal likelihood (grid search)
 *     3. Analytical derivative posterior → μ'(t*), σ'(t*)
 *     4. q_i = sign_i · μ'(t*) / Xv(t*)  with 95% CI = ±1.96 σ'(t*)/Xv
 *
 * The paper (§2.2) uses nested-sampling B-splines for the full posterior;
 * the GP method is the closest practical Bayesian equivalent in-browser.
 */

import type { NoisyMeasurement } from "./simulator";
import { gpEstimateRate } from "./gp";

// ── Types ─────────────────────────────────────────────────────────────────────

export interface MetRaCNoiseConfig {
  Xv_cv:    number;  // VCD coefficient of variation
  Glc_abs:  number;  // σ Glucose  [mM]
  Lac_abs:  number;  // σ Lactate  [mM]
  Gln_abs:  number;  // σ Glutamine[mM]
  Glu_abs:  number;  // σ Glutamate[mM]
  NH4_abs:  number;  // σ Ammonium [mM]
  Tit_cv?:  number;  // VCV Titer  (relative, 0 = not measured)
}

export const DEFAULT_METRAC_NOISE: MetRaCNoiseConfig = {
  Xv_cv:   0.05,
  Glc_abs: 0.30,
  Lac_abs: 0.40,
  Gln_abs: 0.15,
  Glu_abs: 0.08,
  NH4_abs: 0.10,
  Tit_cv:  0.08,
};

export interface RawRateEstimate {
  t: number;
  nearBolus: boolean;   // true → interval spans a feed bolus (rate is unreliable)
  q_Glc: number; q_Glc_sd: number;
  q_Lac: number; q_Lac_sd: number;
  q_Gln: number; q_Gln_sd: number;
  q_Glu: number; q_Glu_sd: number;
  q_NH4: number; q_NH4_sd: number;
  q_p:   number; q_p_sd:   number;  // product-specific rate (0 if Tit not measured)
}

export interface SmoothedRate {
  t: number;
  q_Glc: number; q_Glc_lo95: number; q_Glc_hi95: number;
  q_Lac: number; q_Lac_lo95: number; q_Lac_hi95: number;
  q_Gln: number; q_Gln_lo95: number; q_Gln_hi95: number;
  q_Glu: number; q_Glu_lo95: number; q_Glu_hi95: number;
  q_NH4: number; q_NH4_lo95: number; q_NH4_hi95: number;
  q_p:   number; q_p_lo95:   number; q_p_hi95:   number;
}

// ── Bolus detection ────────────────────────────────────────────────────────────

const BOLUS_SIGMA_FACTOR = 10;  // multiply σ by this for near-bolus intervals

function isNearBolus(t_lo: number, t_hi: number, bolusDay: number[]): boolean {
  return bolusDay.some((d) => d > t_lo - 0.05 && d < t_hi + 0.05);
}

// ── Step 1: Finite-difference specific rates ───────────────────────────────────

/**
 * Estimate raw specific rates from noisy concentration measurements.
 *
 * For each point k:
 *   dC/dt ≈ (C_{k+1} − C_{k−1}) / (t_{k+1} − t_{k−1})   [centered]
 *   q = −dC/dt / Xv      (positive = uptake for Glc/Gln/Glu/NH4; Lac sign = net)
 *
 * Intervals that span a feed bolus are flagged and assigned 10× higher σ so the
 * kernel smoother naturally down-weights them.
 *
 * @param measurements  Noisy concentration time-series
 * @param noise         Noise configuration (must match generation noise)
 * @param bolusDay      Feed bolus days  [optional, default = none]
 */
export function estimateRawRates(
  measurements: NoisyMeasurement[],
  noise: MetRaCNoiseConfig,
  bolusDay: number[] = [],
): RawRateEstimate[] {
  const n = measurements.length;
  if (n < 2) return [];

  const results: RawRateEstimate[] = [];

  for (let k = 0; k < n; k++) {
    const prev = measurements[Math.max(0, k - 1)];
    const curr = measurements[k];
    const next = measurements[Math.min(n - 1, k + 1)];

    // Time span and number of noisy points in numerator
    let dt: number, nDiff: number;
    if (k === 0) {
      dt = next.t - curr.t; nDiff = 1;
    } else if (k === n - 1) {
      dt = curr.t - prev.t; nDiff = 1;
    } else {
      dt = next.t - prev.t; nDiff = 2;
    }
    if (dt < 1e-9) continue;

    const nearBolus = isNearBolus(prev.t, next.t, bolusDay);
    const bolusFactor = nearBolus ? BOLUS_SIGMA_FACTOR : 1;

    // Finite differences
    const dGlc = (next.Glc - prev.Glc) / dt;
    const dLac = (next.Lac - prev.Lac) / dt;
    const dGln = (next.Gln - prev.Gln) / dt;
    const dGlu = (next.Glu - prev.Glu) / dt;
    const dNH4 = (next.NH4 - prev.NH4) / dt;

    // Titer rate (if available)
    const hasTit = curr.Tit !== undefined && next.Tit !== undefined && prev.Tit !== undefined;
    const dTit = hasTit
      ? ((next.Tit ?? 0) - (prev.Tit ?? 0)) / dt
      : 0;

    const Xv  = curr.Xv;
    const eps = 1e-4;
    const safeXv = Math.max(Xv, eps);

    const q_Glc = -dGlc / safeXv;
    const q_Lac =  dLac / safeXv;
    const q_Gln = -dGln / safeXv;
    const q_Glu = -dGlu / safeXv;
    const q_NH4 =  dNH4 / safeXv;
    const q_p   = hasTit ?  dTit / safeXv : 0;

    const varXv = (safeXv * noise.Xv_cv) ** 2;

    function qSD(sigma_C: number, q_i: number): number {
      const varDdt = nDiff * sigma_C ** 2 / dt ** 2;
      const varQ   = varDdt / safeXv ** 2 + q_i ** 2 * varXv / safeXv ** 2;
      return bolusFactor * Math.sqrt(Math.max(varQ, 1e-12));
    }

    const titSD = hasTit && noise.Tit_cv
      ? qSD(Math.max(curr.Tit ?? 1, 1) * noise.Tit_cv, q_p)
      : 1e6;   // unmeasured → infinite uncertainty

    results.push({
      t: curr.t,
      nearBolus,
      q_Glc, q_Glc_sd: qSD(noise.Glc_abs, q_Glc),
      q_Lac, q_Lac_sd: qSD(noise.Lac_abs, q_Lac),
      q_Gln, q_Gln_sd: qSD(noise.Gln_abs, q_Gln),
      q_Glu, q_Glu_sd: qSD(noise.Glu_abs, q_Glu),
      q_NH4, q_NH4_sd: qSD(noise.NH4_abs, q_NH4),
      q_p,   q_p_sd:   titSD,
    });
  }

  return results;
}

// ── Step 2: Nadaraya-Watson kernel smoother ────────────────────────────────────

/**
 * Gaussian kernel smoother with inverse-variance weighting.
 *
 *   q̂(t*) = Σ w_k q_k / Σ w_k     where  w_k = K(t*,t_k) / σ²_k
 *   Var(q̂) = 1 / Σ w_k             (posterior variance, conjugate Gaussian)
 */
function kernelSmooth(
  outputTimes: number[],
  rawEstimates: RawRateEstimate[],
  key: keyof RawRateEstimate,
  sdKey: keyof RawRateEstimate,
  h = 1.5,
): { mean: number[]; lo95: number[]; hi95: number[] } {
  return {
    mean:  new Array(outputTimes.length),
    lo95:  new Array(outputTimes.length),
    hi95:  new Array(outputTimes.length),
  };
  // (filled below in-place to avoid allocating three times)
  // ↑ dummy; replaced with the real implementation:
}

function kernelSmoothImpl(
  outputTimes: number[],
  rawEstimates: RawRateEstimate[],
  key: keyof RawRateEstimate,
  sdKey: keyof RawRateEstimate,
  h = 1.5,
): { mean: number[]; lo95: number[]; hi95: number[] } {
  const mean: number[] = [];
  const lo95: number[] = [];
  const hi95: number[] = [];

  for (const t_star of outputTimes) {
    let wSum = 0, wySum = 0, w2s2Sum = 0;
    for (const r of rawEstimates) {
      const K  = Math.exp(-((t_star - r.t) ** 2) / (2 * h * h));
      const sd = Math.max(r[sdKey] as number, 1e-12);
      const w  = K / (sd * sd);
      wSum    += w;
      wySum   += w * (r[key] as number);
      w2s2Sum += w * w * sd * sd;
    }
    if (wSum < 1e-12) {
      mean.push(0); lo95.push(0); hi95.push(0);
    } else {
      const m   = wySum / wSum;
      const s95 = 1.96 * Math.sqrt(w2s2Sum / (wSum * wSum));
      mean.push(m);
      lo95.push(m - s95);
      hi95.push(m + s95);
    }
  }
  return { mean, lo95, hi95 };
}

// Suppress dead-code warning for the dummy stub
void kernelSmooth;

// ── Step 3: Full MetRaC pipeline ───────────────────────────────────────────────

/**
 * Run the full MetRaC pipeline on noisy measurements.
 *
 * @param measurements  Noisy concentration time-series
 * @param noise         Noise configuration
 * @param outputTimes   Dense time grid for smoothed output
 * @param bandwidth     Kernel bandwidth [days]  (default 1.5, only for "kernel" method)
 * @param bolusDay      Feed bolus days (for discontinuity correction)
 * @param method        "kernel" (default) or "gp" (GP regression with SE kernel)
 */
export function runMetRaC(
  measurements: NoisyMeasurement[],
  noise: MetRaCNoiseConfig,
  outputTimes: number[],
  bandwidth = 1.5,
  bolusDay: number[] = [],
  method: "kernel" | "gp" = "kernel",
): SmoothedRate[] {
  if (method === "gp") {
    return runGPMetRaC(measurements, noise, outputTimes);
  }

  const raw = estimateRawRates(measurements, noise, bolusDay);
  if (raw.length < 2) return [];

  const smooth = (k: keyof RawRateEstimate, sk: keyof RawRateEstimate) =>
    kernelSmoothImpl(outputTimes, raw, k, sk, bandwidth);

  const glc = smooth("q_Glc", "q_Glc_sd");
  const lac = smooth("q_Lac", "q_Lac_sd");
  const gln = smooth("q_Gln", "q_Gln_sd");
  const glu = smooth("q_Glu", "q_Glu_sd");
  const nh4 = smooth("q_NH4", "q_NH4_sd");
  const qp  = smooth("q_p",   "q_p_sd");

  return outputTimes.map((t, i) => ({
    t,
    q_Glc: glc.mean[i], q_Glc_lo95: glc.lo95[i], q_Glc_hi95: glc.hi95[i],
    q_Lac: lac.mean[i], q_Lac_lo95: lac.lo95[i], q_Lac_hi95: lac.hi95[i],
    q_Gln: gln.mean[i], q_Gln_lo95: gln.lo95[i], q_Gln_hi95: gln.hi95[i],
    q_Glu: glu.mean[i], q_Glu_lo95: glu.lo95[i], q_Glu_hi95: glu.hi95[i],
    q_NH4: nh4.mean[i], q_NH4_lo95: nh4.lo95[i], q_NH4_hi95: nh4.hi95[i],
    q_p:   qp.mean[i],  q_p_lo95:   qp.lo95[i],  q_p_hi95:   qp.hi95[i],
  }));
}

// ── GP-based MetRaC pipeline ───────────────────────────────────────────────────

/**
 * GP implementation: fit SE-kernel GP to raw concentration measurements,
 * then compute the analytical derivative posterior to get q values and CIs.
 */
function runGPMetRaC(
  measurements: NoisyMeasurement[],
  noise: MetRaCNoiseConfig,
  outputTimes: number[],
): SmoothedRate[] {
  if (measurements.length < 3) return [];

  const times    = measurements.map((m) => m.t);
  const xvValues = measurements.map((m) => m.Xv);

  const est = (
    values: number[],
    sn: number,
    sign: number,
  ) => gpEstimateRate(times, values, sn, outputTimes, times, xvValues, sign);

  const glc = est(measurements.map((m) => m.Glc), noise.Glc_abs, -1);
  const lac = est(measurements.map((m) => m.Lac), noise.Lac_abs, +1);
  const gln = est(measurements.map((m) => m.Gln), noise.Gln_abs, -1);
  const glu = est(measurements.map((m) => m.Glu), noise.Glu_abs, -1);
  const nh4 = est(measurements.map((m) => m.NH4), noise.NH4_abs, +1);

  // Product-specific rate (titer) — optional
  const hasTit = measurements.some((m) => m.Tit !== undefined) && (noise.Tit_cv ?? 0) > 0;
  const ZERO = new Array<number>(outputTimes.length).fill(0);
  let qp = { mean: ZERO, lo95: ZERO, hi95: ZERO, learnedL: 1.5 };
  if (hasTit) {
    const titValues  = measurements.map((m) => m.Tit ?? 0);
    const meanTit    = titValues.reduce((s, v) => s + v, 0) / titValues.length;
    const titSn      = Math.max(1, meanTit) * (noise.Tit_cv ?? 0.08);
    qp = gpEstimateRate(times, titValues, titSn, outputTimes, times, xvValues, +1);
  }

  return outputTimes.map((t, i) => ({
    t,
    q_Glc: glc.mean[i], q_Glc_lo95: glc.lo95[i], q_Glc_hi95: glc.hi95[i],
    q_Lac: lac.mean[i], q_Lac_lo95: lac.lo95[i], q_Lac_hi95: lac.hi95[i],
    q_Gln: gln.mean[i], q_Gln_lo95: gln.lo95[i], q_Gln_hi95: gln.hi95[i],
    q_Glu: glu.mean[i], q_Glu_lo95: glu.lo95[i], q_Glu_hi95: glu.hi95[i],
    q_NH4: nh4.mean[i], q_NH4_lo95: nh4.lo95[i], q_NH4_hi95: nh4.hi95[i],
    q_p:   qp.mean[i],  q_p_lo95:   qp.lo95[i],  q_p_hi95:   qp.hi95[i],
  }));
}
