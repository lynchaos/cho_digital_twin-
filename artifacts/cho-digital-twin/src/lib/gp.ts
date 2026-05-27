/**
 * Gaussian Process regression for MetRaC rate estimation.
 *
 * Fits a GP with squared-exponential (SE/RBF) kernel to noisy concentration
 * measurements, then returns the analytical derivative posterior (dC/dt).
 * Dividing by Xv gives specific metabolic rates with proper Bayesian CIs.
 *
 * Key equations for SE kernel  k(t,t') = sf² exp(−(t−t')²/2l²):
 *   GP posterior mean:          μ(t*)  = K(t*,T) Kₙ⁻¹ y
 *   Derivative posterior mean:  μ'(t*) = ∂K(t*,T)/∂t* · Kₙ⁻¹ y
 *   Derivative posterior var:   v'(t*) = sf²/l² − ∂K(t*,T)/∂t* · Kₙ⁻¹ · ∂K(T,t*)/∂t*
 *
 * Hyperparameter l (length-scale) is optimised by grid-search on the log
 * marginal likelihood.  sf is estimated from the data range; sn is supplied
 * by the caller (known measurement noise σ).
 */

// ── Linear algebra helpers ────────────────────────────────────────────────────

/** In-place Cholesky decomposition: K = L Lᵀ (K must be positive-definite). */
function cholesky(K: number[][]): number[][] {
  const n = K.length;
  const L: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = K[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      L[i][j] = j < i
        ? s / Math.max(L[j][j], 1e-12)
        : Math.sqrt(Math.max(s, 1e-10));
    }
  }
  return L;
}

/** Forward substitution: solve L x = b. */
function fwdSub(L: number[][], b: number[]): number[] {
  const n = b.length;
  const x = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = b[i];
    for (let k = 0; k < i; k++) s -= L[i][k] * x[k];
    x[i] = s / L[i][i];
  }
  return x;
}

/** Back substitution: solve Lᵀ x = b. */
function bwdSub(L: number[][], b: number[]): number[] {
  const n = b.length;
  const x = new Array<number>(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k];
    x[i] = s / L[i][i];
  }
  return x;
}

/** Solve K x = b where K = L Lᵀ (Cholesky). */
function cholSolve(L: number[][], b: number[]): number[] {
  return bwdSub(L, fwdSub(L, b));
}

function dot(a: number[], b: number[]): number {
  return a.reduce((s, x, i) => s + x * b[i], 0);
}

// ── SE kernel and its derivative ──────────────────────────────────────────────

function seCov(t1: number, t2: number, l: number, sf: number): number {
  const r = (t1 - t2) / l;
  return sf * sf * Math.exp(-0.5 * r * r);
}

/** ∂k(t*,t)/∂t*  — derivative of SE kernel w.r.t. first argument. */
function seCovDeriv(tStar: number, t: number, l: number, sf: number): number {
  return -(tStar - t) / (l * l) * seCov(tStar, t, l, sf);
}

// ── Kernel matrix builders ────────────────────────────────────────────────────

function buildKTT(times: number[], l: number, sf: number, sn: number): number[][] {
  const n = times.length;
  return Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (__, j) =>
      seCov(times[i], times[j], l, sf) + (i === j ? sn * sn + 1e-6 : 0)
    )
  );
}

// ── Log marginal likelihood ───────────────────────────────────────────────────

function logML(L: number[][], alpha: number[], y: number[]): number {
  let logDet = 0;
  for (let i = 0; i < L.length; i++) logDet += Math.log(Math.max(L[i][i], 1e-12));
  return -0.5 * dot(y, alpha) - logDet - 0.5 * y.length * Math.log(2 * Math.PI);
}

// ── Hyperparameter optimisation ───────────────────────────────────────────────

const L_GRID = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0];

function optimiseL(times: number[], values: number[], sf: number, sn: number): number {
  let bestL = 1.5, bestScore = -Infinity;
  for (const l of L_GRID) {
    try {
      const K = buildKTT(times, l, sf, sn);
      const Lc = cholesky(K);
      const alpha = cholSolve(Lc, values);
      const score = logML(Lc, alpha, values);
      if (Number.isFinite(score) && score > bestScore) {
        bestScore = score;
        bestL = l;
      }
    } catch { /* ill-conditioned — skip */ }
  }
  return bestL;
}

// ── GP model ──────────────────────────────────────────────────────────────────

interface GPModel {
  times: number[];
  alpha: number[];  // K⁻¹ y
  Lc:    number[][];  // Cholesky factor (reused for variance)
  l: number; sf: number; sn: number;
}

function fitGP(times: number[], values: number[], sn: number): GPModel {
  const lo = Math.min(...values), hi = Math.max(...values);
  const sf = Math.max(0.01, (hi - lo) * 0.5);
  const l  = optimiseL(times, values, sf, sn);
  const K  = buildKTT(times, l, sf, sn);
  const Lc = cholesky(K);
  const alpha = cholSolve(Lc, values);
  return { times, alpha, Lc, l, sf, sn };
}

interface DerivPred { mean: number; std: number; }

function predictDeriv(model: GPModel, tStar: number): DerivPred {
  const { times, alpha, Lc, l, sf } = model;

  // Derivative cross-covariance: dks_i = ∂k(t*,t_i)/∂t*
  const dks = times.map(t => seCovDeriv(tStar, t, l, sf));

  // Derivative posterior mean: μ'(t*) = dks @ α
  const dmean = dot(dks, alpha);

  // Derivative posterior variance: v'(t*) = sf²/l² − dks @ K⁻¹ @ dksᵀ
  const Kinv_dks = cholSolve(Lc, dks);
  const priorVar = sf * sf / (l * l);
  const postVar  = Math.max(0, priorVar - dot(dks, Kinv_dks));

  return { mean: dmean, std: Math.sqrt(postVar) };
}

// ── Linear interpolation helper ───────────────────────────────────────────────

function lerp(xs: number[], ys: number[], x: number): number {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[xs.length - 1];
  const i = xs.findIndex(v => v >= x) - 1;
  const t = (x - xs[i]) / (xs[i + 1] - xs[i]);
  return ys[i] + t * (ys[i + 1] - ys[i]);
}

// ── Public API ────────────────────────────────────────────────────────────────

export interface GPRateResult {
  mean:    number[];
  lo95:    number[];
  hi95:    number[];
  learnedL: number;  // optimised length-scale [days], for display
}

/**
 * Estimate a specific metabolic rate via GP derivative posterior.
 *
 * @param measTimes   Measurement time points [days]
 * @param measValues  Noisy concentration measurements [mM or mg/L]
 * @param sn          Measurement noise σ (known from instrument spec)
 * @param outputTimes Dense time grid for output
 * @param xvTimes     Time points at which Xv is known (measurement times)
 * @param xvValues    Viable cell density at xvTimes [Mc/mL]
 * @param sign        +1 secretion (Lac, NH4, Tit), −1 uptake (Glc, Gln, Glu)
 */
export function gpEstimateRate(
  measTimes:   number[],
  measValues:  number[],
  sn:          number,
  outputTimes: number[],
  xvTimes:     number[],
  xvValues:    number[],
  sign:        number,
): GPRateResult {
  const EMPTY: GPRateResult = {
    mean:    new Array<number>(outputTimes.length).fill(0),
    lo95:    new Array<number>(outputTimes.length).fill(0),
    hi95:    new Array<number>(outputTimes.length).fill(0),
    learnedL: 1.5,
  };
  if (measTimes.length < 3) return EMPTY;

  const model = fitGP(measTimes, measValues, sn);

  const mean: number[] = [];
  const lo95: number[] = [];
  const hi95: number[] = [];

  for (const t of outputTimes) {
    const pred = predictDeriv(model, t);
    const xv   = Math.max(lerp(xvTimes, xvValues, t), 1e-4);

    const meanQ = sign * pred.mean / xv;
    const stdQ  = pred.std / xv;

    mean.push(meanQ);
    lo95.push(meanQ - 1.96 * stdQ);
    hi95.push(meanQ + 1.96 * stdQ);
  }

  return { mean, lo95, hi95, learnedL: model.l };
}
