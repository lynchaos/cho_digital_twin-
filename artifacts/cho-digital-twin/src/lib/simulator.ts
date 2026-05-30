
/**
 * CHO Cell Culture Fed-Batch Simulator
 *
 * Integrates the coupled ODE system (biomass + FLEX) with:
 *  • Nutrient-coupled μ_net (substitute for paper's §2.3 NN)
 *  • Luedeking-Piret product model (§2.4.2 extension)
 *  • Bolus feeding with exact volume mixing
 *  • Optional noisy "virtual measurements" for MetRaC input generation
 */

import { rk4Step } from "./ode-solver";
import { choODE, type FeedInput } from "./models";
import {
  sigmaBaseline, nutrientCoupledMuNet, nnMuNet, type SigmoidComponent, type MuNetMode,
} from "./growth-rate";
import type { NNWeights } from "./neural-net";
import {
  DEFAULT_MODEL_PARAMS, DEFAULT_NUTRIENT_COUPLING,
  type ModelParams, type NutrientCouplingParams,
} from "./parameters";

// ── Types ────────────────────────────────────────────────────────────────────

export interface FeedBolus {
  day: number;
  volumeFraction: number;
  Glc_feed: number; Gln_feed: number; Glu_feed: number;
}

export interface InitialConditions {
  Xv: number; Xd: number; Xl: number; B: number;
  Glc: number; Lac: number; Gln: number; Glu: number; NH4: number; Tit: number;
}

export interface SimulationConfig {
  initialConditions: InitialConditions;
  initialVolume: number;           // [mL]
  feedBoluses: FeedBolus[];
  muNetComponents: SigmoidComponent[];
  nutrientCoupling: NutrientCouplingParams;
  modelParams: ModelParams;
  runDays: number;
  outputInterval: number;          // [days]
  muNetMode?: MuNetMode;           // default "nutrient-coupled"
  nnWeights?: NNWeights | null;    // required when muNetMode === "surrogate-nn"
}

export interface TimePoint {
  t: number; Xv: number; Xd: number; Xl: number; B: number;
  Glc: number; Lac: number; Gln: number; Glu: number; NH4: number; Tit: number;
  mu_net: number; mu_eff: number; kd: number; kl: number; volume: number;
  // specific rates (for MetRaC comparison)
  q_Glc: number; q_Lac: number; q_Gln: number; q_Glu: number; q_NH4: number;
  q_p: number;   // product-specific rate [mg·L⁻¹·(Mc/mL)⁻¹·day⁻¹]
}

export interface NoisyMeasurement {
  t: number;
  Xv: number; Glc: number; Lac: number; Gln: number; Glu: number; NH4: number;
  Tit?: number;    // product titer [mg/L] — optional (if Tit_cv > 0)
  volume: number;
  // true (noiseless) values for comparison
  Xv_true: number; Glc_true: number; Lac_true: number;
  Gln_true: number; Glu_true: number; NH4_true: number;
  Tit_true?: number;
}

// ── Defaults ─────────────────────────────────────────────────────────────────

export const DEFAULT_INITIAL_CONDITIONS: InitialConditions = {
  Xv: 0.30, Xd: 0, Xl: 0, B: 0,
  Glc: 26.0, Lac: 0.5, Gln: 8.0, Glu: 2.0, NH4: 0.2, Tit: 0,
};

export const DEFAULT_FEED_BOLUSES: FeedBolus[] = [
  { day:  3, volumeFraction: 0.05, Glc_feed: 100, Gln_feed: 40, Glu_feed: 8 },
  { day:  5, volumeFraction: 0.05, Glc_feed: 100, Gln_feed: 40, Glu_feed: 8 },
  { day:  7, volumeFraction: 0.05, Glc_feed: 100, Gln_feed: 40, Glu_feed: 8 },
  { day:  9, volumeFraction: 0.05, Glc_feed: 100, Gln_feed: 40, Glu_feed: 8 },
  { day: 11, volumeFraction: 0.05, Glc_feed: 100, Gln_feed: 40, Glu_feed: 8 },
];

// ── Internal helpers ──────────────────────────────────────────────────────────

function applyBolus(y: number[], volume: number, bolus: FeedBolus): [number[], number] {
  const dV = bolus.volumeFraction * volume;
  const V_new = volume + dV;
  const f = volume / V_new;
  const yn = [...y];
  yn[0] = y[0]*f; yn[1] = y[1]*f; yn[2] = y[2]*f; yn[3] = y[3]*f;
  yn[4] = (y[4]*volume + bolus.Glc_feed*dV) / V_new;
  yn[5] = y[5]*f;
  yn[6] = (y[6]*volume + bolus.Gln_feed*dV) / V_new;
  yn[7] = (y[7]*volume + bolus.Glu_feed*dV) / V_new;
  yn[8] = y[8]*f;
  yn[9] = y[9]*f;
  return [yn, V_new];
}

function computeMuNet(
  t: number, y: number[],
  muComponents: SigmoidComponent[],
  nc: NutrientCouplingParams,
  mode: MuNetMode = "nutrient-coupled",
  nnWeights: NNWeights | null | undefined = null,
): number {
  const [, , , , Glc, Lac, Gln, , NH4] = y;
  if (mode === "surrogate-nn") {
    return nnMuNet(t, Glc, Gln, Lac, NH4, nc, muComponents, nnWeights ?? null);
  }
  const mu_base = sigmaBaseline(t, muComponents);
  if (mode === "sigmoid") return mu_base;
  return nutrientCoupledMuNet(mu_base, Glc, Gln, Lac, NH4, nc);
}

// ── Simulator ─────────────────────────────────────────────────────────────────

import {
  qGlc as modelQGlc, qLacNet as modelQLac,
  qGln as modelQGln, qGlu as modelQGlu, qNH4 as modelQNH4,
  deathRate, lysisRate,
} from "./models";

export function runSimulation(config: SimulationConfig): TimePoint[] {
  const { initialConditions: ic, initialVolume, feedBoluses,
          muNetComponents, nutrientCoupling: nc, modelParams: mp,
          runDays, outputInterval,
          muNetMode = "nutrient-coupled", nnWeights = null } = config;

  const dt = 0.005;
  const nOut = Math.round(runDays / outputInterval) + 1;
  const outputTimes = Array.from({ length: nOut }, (_, i) => i * outputInterval);
  const sortedBoluses = [...feedBoluses].sort((a, b) => a.day - b.day);

  let y: number[] = [ic.Xv, ic.Xd, ic.Xl, ic.B, ic.Glc, ic.Lac, ic.Gln, ic.Glu, ic.NH4, ic.Tit];
  let volume = initialVolume;
  const results: TimePoint[] = [];
  let t = 0, bolusIdx = 0;

  const record = (t: number, y: number[], vol: number) => {
    const mu_net = computeMuNet(t, y, muNetComponents, nc, muNetMode, nnWeights);
    const kd_v   = deathRate(y[3], mp);
    const kl_v   = lysisRate(y[2], mp);
    const mu_eff = mu_net + kd_v;
    const q_glc  = modelQGlc(mu_eff, y[4], y[5], mp);
    const q_lac  = modelQLac(q_glc, y[5], mp);
    const q_glu  = modelQGlu(mu_eff, y[7], mp);
    const q_gln  = modelQGln(y[6], mp);
    const q_nh4  = modelQNH4(q_glu, q_gln, mp);
    results.push({
      t,
      Xv: Math.max(0, y[0]),  Xd: Math.max(0, y[1]),  Xl: Math.max(0, y[2]),
      B:  Math.max(0, y[3]),  Glc: Math.max(0, y[4]), Lac: Math.max(0, y[5]),
      Gln: Math.max(0, y[6]), Glu: Math.max(0, y[7]), NH4: Math.max(0, y[8]),
      Tit: Math.max(0, y[9]),
      mu_net, mu_eff, kd: kd_v, kl: kl_v, volume: vol,
      q_Glc: q_glc, q_Lac: q_lac, q_Gln: q_gln, q_Glu: q_glu, q_NH4: q_nh4,
      q_p: mp.q_p_growth * Math.max(0, mu_net) + mp.q_p,
    });
  };

  record(0, y, volume);

  for (let oi = 1; oi < outputTimes.length; oi++) {
    const t_out = outputTimes[oi];

    while (t < t_out - 1e-9) {
      while (bolusIdx < sortedBoluses.length && sortedBoluses[bolusIdx].day <= t + 1e-9) {
        [y, volume] = applyBolus(y, volume, sortedBoluses[bolusIdx++]);
      }
      const stepDt = Math.min(dt, t_out - t);
      const feed: FeedInput = { F: 0, V: volume, Glc_feed: 0, Gln_feed: 0, Glu_feed: 0 };
      y = rk4Step(
        (tt, yy) => choODE(tt, yy, computeMuNet(tt, yy, muNetComponents, nc, muNetMode, nnWeights), feed, mp),
        t, y, stepDt,
      );
      t += stepDt;
    }
    while (bolusIdx < sortedBoluses.length && sortedBoluses[bolusIdx].day <= t + 1e-9) {
      [y, volume] = applyBolus(y, volume, sortedBoluses[bolusIdx++]);
    }
    record(t_out, y, volume);
  }
  return results;
}

export function defaultConfig(overrides?: Partial<SimulationConfig>): SimulationConfig {
  return {
    initialConditions: DEFAULT_INITIAL_CONDITIONS,
    initialVolume: 14,
    feedBoluses: DEFAULT_FEED_BOLUSES,
    muNetComponents: [
      { a:  0.85, b: 1.4, c:  1.5 },
      { a: -0.95, b: 0.9, c:  7.0 },
      { a:  0.25, b: 0.5, c:  4.0 },
      { a: -0.18, b: 1.2, c: 11.0 },
    ],
    nutrientCoupling: { ...DEFAULT_NUTRIENT_COUPLING },
    modelParams: { ...DEFAULT_MODEL_PARAMS },
    runDays: 14,
    outputInterval: 0.1,
    ...overrides,
  };
}

// ── Noisy measurement generator (for MetRaC input) ───────────────────────────

/** Box-Muller normal random number generator */
function randn(mean: number, sd: number): number {
  const u1 = Math.random(), u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1 + 1e-10)) * Math.cos(2 * Math.PI * u2);
  return mean + sd * z;
}

export interface MeasurementNoiseConfig {
  Xv_cv: number;   // coefficient of variation for VCD (e.g. 0.05 = 5%)
  Glc_abs: number; // absolute noise σ for Glc [mM]
  Lac_abs: number; // absolute noise σ for Lac [mM]
  Gln_abs: number; // absolute noise σ for Gln [mM]
  Glu_abs: number; // absolute noise σ for Glu [mM]
  NH4_abs: number; // absolute noise σ for NH4 [mM]
  sampleEvery: number; // sampling interval [days]
  Tit_cv?: number; // relative noise for Titer (CV, e.g. 0.08 = 8%); omit = not measured
}

export const DEFAULT_NOISE_CONFIG: MeasurementNoiseConfig = {
  Xv_cv: 0.05,
  Glc_abs: 0.3,
  Lac_abs: 0.4,
  Gln_abs: 0.15,
  Glu_abs: 0.08,
  NH4_abs: 0.1,
  sampleEvery: 1.0,
  Tit_cv: 0.08,
};

/**
 * Generate synthetic "experimental" measurements by subsampling the ODE
 * output and adding Gaussian noise.  Used as MetRaC input.
 */
export function generateNoisyMeasurements(
  simResults: TimePoint[],
  noise: MeasurementNoiseConfig,
): NoisyMeasurement[] {
  // Pick output points at the requested sample interval
  const sampled = simResults.filter((r) => {
    const frac = r.t / noise.sampleEvery;
    return Math.abs(frac - Math.round(frac)) < 0.01;
  });

  return sampled.map((r) => ({
    t:        r.t,
    volume:   r.volume,
    Xv_true:  r.Xv,
    Glc_true: r.Glc,
    Lac_true: r.Lac,
    Gln_true: r.Gln,
    Glu_true: r.Glu,
    NH4_true: r.NH4,
    Tit_true: r.Tit,
    Xv:  Math.max(0, randn(r.Xv,  r.Xv  * noise.Xv_cv)),
    Glc: Math.max(0, randn(r.Glc, noise.Glc_abs)),
    Lac: Math.max(0, randn(r.Lac, noise.Lac_abs)),
    Gln: Math.max(0, randn(r.Gln, noise.Gln_abs)),
    Glu: Math.max(0, randn(r.Glu, noise.Glu_abs)),
    NH4: Math.max(0, randn(r.NH4, noise.NH4_abs)),
    Tit: noise.Tit_cv && noise.Tit_cv > 0
      ? Math.max(0, randn(r.Tit, Math.max(r.Tit * noise.Tit_cv, 5)))
      : undefined,
  }));
}

/** Export TimePoint[] to CSV string */
export function exportToCsv(results: TimePoint[]): string {
  const header = [
    "t_day","Xv_Mc_per_mL","Xd_Mc_per_mL","Xl_Mc_per_mL","B",
    "Glc_mM","Lac_mM","Gln_mM","Glu_mM","NH4_mM","Tit_mg_per_L",
    "mu_net_per_day","mu_eff_per_day","kd_per_day","kl_per_day",
    "volume_mL","q_Glc","q_Lac","q_Gln","q_Glu","q_NH4","q_p",
  ].join(",");
  const rows = results.map((r) =>
    [r.t, r.Xv, r.Xd, r.Xl, r.B, r.Glc, r.Lac, r.Gln, r.Glu, r.NH4, r.Tit,
     r.mu_net, r.mu_eff, r.kd, r.kl, r.volume,
     r.q_Glc, r.q_Lac, r.q_Gln, r.q_Glu, r.q_NH4, r.q_p]
    .map((v) => v.toFixed(6)).join(","),
  );
  return [header, ...rows].join("\n");
}
