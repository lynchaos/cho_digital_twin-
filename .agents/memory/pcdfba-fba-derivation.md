---
name: PC-dFBA analytical FBA derivation
description: Sign conventions and analytical solution for the condensed CHO 10-metabolite / 16-reaction FBA network in cho-network.ts / fba-solver.ts
---

## Network

Files: `artifacts/cho-digital-twin/src/lib/cho-network.ts` (S matrix + SVG positions),
`artifacts/cho-digital-twin/src/lib/fba-solver.ts` (analytical solver + PCA),
`artifacts/cho-digital-twin/src/pages/PcdFBAPage.tsx` (UI).

10 intracellular metabolites: PYR, LAC, AcCoA, OAA, CIT, AKG, SUCC, MAL, GLU, GLN
4 exchange + 12 internal = 16 reactions total.

## Sign Conventions (critical — match TimePoint in simulator.ts)

| Variable | Positive means |
|---|---|
| EX_GLC = q_Glc | cells consume glucose (uptake) |
| EX_LAC = q_Lac | cells produce lactate (secretion into broth) |
| EX_GLN = q_Gln | cells consume glutamine (uptake) |
| EX_GLU = q_Glu | cells consume glutamate from broth (UPTAKE, NOT secretion) |

**EX_GLU is uptake-positive.** The S matrix column for EX_GLU has GLU = +1 (uptake increases intracellular pool). This was deliberately set to match `q_Glu` from `TimePoint` which records consumption rate. Note: transamination-driven Glu secretion (typical in early CHO culture) is NOT captured by this condensed model.

## Analytical Solution (parsimonious FBA)

Derived from S·v = 0 with parsimony choice ME = GLUD (minimise anaplerosis → PCX = 0):

```
GLS  = EX_GLN
GLUD = EX_GLN + EX_GLU       (total Glu catabolism = Gln-derived + direct uptake)
PDH  = 2·EX_GLC − EX_LAC + GLUD   (fully determined by mass balance)
CS   = ICDH = PDH
AKGDH = SDH_FH = PDH + GLUD
ME   = GLUD   (parsimony)
MDH  = PDH    (follows from ME = GLUD)
PCX  = 0      (follows from ME = GLUD)
LDH_f = max(0, EX_LAC)
LDH_r = max(0, -EX_LAC)
```

**Why:** The 10-metabolite mass balance leaves one free parameter (ME vs PCX split). Parsimony assigns ME = minimum needed (= GLUD to keep PCX ≥ 0). This gives a unique analytical solution with zero LP overhead.

**Feasibility flag:** PDH < 0 → infeasible (clamped to 0 with warning). Occurs if lactate re-consumption flux exceeds 2·EX_GLC + GLUD.

## PCA

`computeFluxPCA(FBAResult[])` in fba-solver.ts. Uses 12 internal reactions (excludes exchange). Power iteration for PC1, deflated power iteration for PC2. Scores represent the culture-trajectory in metabolic state space.
