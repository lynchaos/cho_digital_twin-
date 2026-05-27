/**
 * Analytical FBA solver for the condensed CHO central carbon metabolism network.
 *
 * After fixing the 4 measured exchange rates, the 10-metabolite mass balance
 * uniquely determines most internal fluxes. The remaining degree of freedom
 * (malic enzyme ME vs. pyruvate carboxylase PCX anaplerosis split) is resolved
 * by the parsimonious choice: ME = GLUD, PCX = 0.
 *
 * Derived analytically from S·v = 0:
 *   GLS  = EX_GLN
 *   GLUD = GLS + EX_GLU  (total Glu catabolism = Gln-derived + direct import)
 *   PDH  = 2·EX_GLC − EX_LAC + GLUD   (mass balance, fully determined)
 *   CS = ICDH = PDH
 *   AKGDH = SDH = PDH + GLUD
 *   ME  = GLUD  (parsimony: minimise anaplerosis)
 *   MDH = PDH   (follows from ME = GLUD)
 *   PCX = 0     (follows from ME = GLUD)
 *   LDH_f = max(0, EX_LAC),  LDH_r = max(0, −EX_LAC)
 *
 * Exchange flux sign conventions (match TimePoint in simulator.ts):
 *   q_Glc > 0  = glucose uptake   q_Lac > 0  = lactate secretion
 *   q_Gln > 0  = glutamine uptake  q_Glu > 0  = glutamate uptake from broth
 */

import type { RxnId } from "./cho-network";

export type FBAFluxes = Record<RxnId, number>;

export interface FBAResult {
  fluxes: FBAFluxes;
  status: "feasible" | "infeasible";
  warnings: string[];
}

export function runFBA(
  q_Glc: number,  // glucose uptake rate   (positive = cells consume Glc)
  q_Lac: number,  // net lactate exchange  (positive = cells produce/secrete Lac)
  q_Gln: number,  // glutamine uptake rate (positive = cells consume Gln)
  q_Glu: number,  // glutamate uptake rate (positive = cells consume Glu from broth)
): FBAResult {
  const warnings: string[] = [];

  const p = Math.max(0, q_Glc);
  const l = q_Lac;                       // net lactate (may be negative = re-consumption)
  const q = Math.max(0, q_Gln);
  const g = Math.max(0, q_Glu);

  const GLS  = q;
  const GLUD = Math.max(0, GLS + g);

  const PDH_raw = 2 * p - l + GLUD;
  const PDH  = Math.max(0, PDH_raw);
  if (PDH_raw < -0.01) warnings.push("PDH < 0: excess lactate exceeds glycolytic supply");

  const CS    = PDH;
  const ICDH  = PDH;
  const AKGDH = PDH + GLUD;
  const SDH   = AKGDH;
  const ME    = GLUD;
  const PCX   = 0;
  const MDH   = Math.max(0, SDH - ME);

  const LDH_f = Math.max(0,  l);
  const LDH_r = Math.max(0, -l);

  const feasible =
    PDH >= -1e-9 && MDH >= -1e-9 && GLS >= -1e-9 && GLUD >= -1e-9;

  return {
    fluxes: {
      EX_GLC: p, EX_LAC: l, EX_GLN: q, EX_GLU: g,
      LDH_f, LDH_r, PDH, CS, ICDH, AKGDH, SDH_FH: SDH, MDH, ME, PCX, GLS, GLUD,
    },
    status: feasible ? "feasible" : "infeasible",
    warnings,
  };
}

// ── PCA helpers ───────────────────────────────────────────────────────────────

function transpose(A: number[][]): number[][] {
  const m = A.length, n = A[0].length;
  const T: number[][] = Array.from({ length: n }, () => new Array<number>(m).fill(0));
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) T[j][i] = A[i][j];
  return T;
}
function dot(a: number[], b: number[]): number {
  return a.reduce((s, x, i) => s + x * b[i], 0);
}
function norm(v: number[]): number { return Math.sqrt(dot(v, v)); }
function normalize(v: number[]): number[] {
  const n = norm(v);
  return n < 1e-14 ? v : v.map(x => x / n);
}

function powerIterate(C: number[][], iters = 300): [number[], number] {
  const n = C.length;
  let v = normalize(Array.from({ length: n }, (_, i) => Math.sin(i + 1)));
  for (let k = 0; k < iters; k++) {
    const w: number[] = new Array<number>(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) w[i] += C[i][j] * v[j];
    const len = norm(w);
    if (len < 1e-14) break;
    v = normalize(w);
  }
  const Av: number[] = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Av[i] += C[i][j] * v[j];
  return [v, dot(v, Av)];
}

// Internal reactions included in PCA (omit exchange fluxes)
const PCA_RXNS: RxnId[] = [
  "PDH","CS","ICDH","AKGDH","SDH_FH","MDH","ME","PCX","GLS","GLUD","LDH_f","LDH_r",
];

export interface PCAResult {
  scores: [number, number][];     // [n_points × 2] PC scores
  loadings: [number[], number[]]; // [PC1_loadings, PC2_loadings]
  explained: [number, number];    // fraction of variance explained
  rxnOrder: string[];
}

export function computeFluxPCA(results: FBAResult[]): PCAResult {
  const EMPTY: PCAResult = {
    scores: [], loadings: [[], []], explained: [0, 0], rxnOrder: PCA_RXNS,
  };
  const n = results.length;
  if (n < 4) return EMPTY;

  const p = PCA_RXNS.length;
  const X: number[][] = results.map(r => PCA_RXNS.map(rxn => r.fluxes[rxn] ?? 0));

  // Centre
  const mean = new Array<number>(p).fill(0);
  for (const row of X) row.forEach((v, j) => { mean[j] += v / n; });
  const Xc = X.map(row => row.map((v, j) => v - mean[j]));

  // Covariance C = Xc.T @ Xc / (n-1)
  const XT = transpose(Xc);
  const C: number[][] = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (__, j) => dot(XT[i], XT[j]) / (n - 1))
  );
  const totalVar = C.reduce((s, row, i) => s + row[i], 0);
  if (totalVar < 1e-14) return EMPTY;

  // PC1
  const [pc1, lam1] = powerIterate(C);

  // Deflate for PC2
  const C2 = C.map((row, i) =>
    row.map((v, j) => v - lam1 * pc1[i] * pc1[j])
  );
  const [pc2, lam2] = powerIterate(C2);

  const scores: [number, number][] = Xc.map(
    row => [dot(row, pc1), dot(row, pc2)] as [number, number]
  );
  return {
    scores,
    loadings: [pc1, pc2],
    explained: [lam1 / totalVar, lam2 / totalVar],
    rxnOrder: PCA_RXNS,
  };
}
