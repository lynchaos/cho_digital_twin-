# CHO Digital Twin

Interactive replication of the Richelle et al. (2025) hybrid CHO cell culture modeling framework (bioRxiv 2025.11.24.690194) as a React web app.

## Run & Operate

- `pnpm --filter @workspace/cho-digital-twin run dev` — run the main app (uses `PORT` env var)
- `pnpm run typecheck` — full typecheck across all packages
- `pnpm run build` — typecheck + build all packages

## Stack

- pnpm workspaces, Node.js 24, TypeScript 5.9
- Frontend: React 19 + Vite 7 + Tailwind CSS + Recharts
- Routing: Wouter
- No backend (pure client-side computation)

## Where things live

- `artifacts/cho-digital-twin/src/` — main app source
  - `pages/` — one file per tab (SimulatorPage, EquationsPage, ParametersPage, MetRaCPage, SweepPage, PcdFBAPage)
  - `lib/models.ts` — all ODE rate equations (Eqs. 1–26 + Gln/Glu transamination)
  - `lib/simulator.ts` — ODE integrator (RK45-like), defaultConfig, runSimulation
  - `lib/cho-network.ts` — condensed 10-met/16-rxn CHO S matrix + SVG node positions
  - `lib/fba-solver.ts` — analytical FBA solver + PCA (power iteration)
  - `lib/growth-rate.ts` — μ_net model switcher (Sigmoid / Monod proxy / Surrogate NN)
  - `lib/neural-net.ts` — surrogate NN (standalone math)
  - `lib/parameters.ts` — Table 1 parameters (all exact from paper)

## Architecture decisions

- All ODE rate equations are exact replications of Eqs. 8–26 from the paper. Table 1 parameters match to published values.
- μ_net has three selectable modes: Sigmoid (exact paper Eq. 14), Monod-proxy (realistic VCD), Surrogate NN (auto-calibrated from Monod proxy).
- PC-dFBA uses an analytical mass-balance solver (no LP library needed) — the condensed 10-metabolite network has only 1 free variable (ME/PCX split), resolved by parsimony.
- EX_GLU sign convention: uptake-positive (matches q_Glu in TimePoint which records consumption rate). See memory/pcdfba-fba-derivation.md.

## Product

Seven tabs:
1. **Simulator** — interactive fed-batch ODE simulation with all controls from the paper
2. **Equations** — rendered LaTeX of all paper equations (Eqs. 1–26 + PC-dFBA Eqs. 27–33)
3. **Parameters** — Table 1 interactive editor with paper default values
4. **MetRaC** — metabolic rate calculator with Bayesian CIs and q_p chart
5. **Sweep** — parameter sensitivity sweep across any model parameter
6. **PC-dFBA** — intracellular flux analysis: SVG flux map + time series + PC trajectory analysis
7. **About** — replication status table

## User preferences

_Populate as you build — explicit user instructions worth remembering across sessions._

## Gotchas

- Do NOT run `pnpm dev` at workspace root — it has no dev script; use the workflow.
- `pnpm run typecheck` (not `build`) for local type-checking of artifacts.
- The surrogate NN needs explicit calibration call in growth-rate.ts; do not import training logic in neural-net.ts (circular import risk).
- EX_GLU in cho-network.ts S matrix has GLU row = +1 (uptake-positive); do NOT change to -1.

## Pointers

- See the `pnpm-workspace` skill for workspace structure, TypeScript setup, and package details
