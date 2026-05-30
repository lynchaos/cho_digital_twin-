# CHO Digital Twin

**An interactive, browser-based digital twin for Chinese Hamster Ovary (CHO) cell culture bioprocessing.**

A live implementation of the hybrid mechanistic + machine learning + genome-scale metabolic modeling framework published by Richelle et al. (2025) and Antonakoudis & Richelle (2026). Run full 14-day fed-batch simulations, explore intracellular flux distributions, and compress a 6,663-reaction genome-scale model — all in your browser, with no installation required.

---

## What It Does

CHO cells are the workhorse of biopharmaceutical manufacturing. Understanding how nutrients, metabolic byproducts, and feeding strategies affect cell growth and monoclonal antibody (mAb) titer is critical — but wet-lab experiments are slow and expensive. This digital twin compresses weeks of experimentation into seconds of simulation by combining:

- **26 coupled ODE equations** describing viable cell density, metabolite concentrations, and product titer in a fed-batch bioreactor
- **Three μ_net (net growth rate) models**: sigmoid baseline, nutrient-coupled Monod proxy, and a surrogate neural network auto-calibrated at load time
- **PC-dFBA** (Principal Component dynamic Flux Balance Analysis) for visualising intracellular metabolic flux distributions
- **MetRaC** (Bayesian metabolic rate calculator) for computing specific exchange rates with credible intervals from concentration time series
- **GEM Reduction Pipeline** — the full 5-step Antonakoudis & Richelle (2026) algorithm that compresses iCHOv1 from 6,663 → 503 reactions (92.5%) in ~22 seconds

---

## Application Tabs

| Tab | Description |
|-----|-------------|
| **Simulator** | 14-day fed-batch ODE simulation with interactive sliders for initial conditions, growth parameters, nutrient coupling constants, feed strategy, and μ_net mode. Six live charts + KPI summary. Runs in <50 ms. |
| **Equations** | All 33 model equations (Eqs. 1–26 kinetic core, Eqs. 27–33 PC-dFBA) rendered in LaTeX notation, grouped by biological subsystem with paper cross-references. |
| **Parameters** | Interactive editor for all 40+ Table 1 parameters (symbol, value, units, equation number). Changes propagate to the Simulator tab. One-click reset to paper defaults. |
| **MetRaC** | Gaussian Process-based metabolic rate calculator. Computes specific exchange rates (q_Glc, q_Lac, q_Gln, q_Glu, q_NH₄, q_mAb) with 95% Bayesian credible intervals from user-uploaded concentration profiles. |
| **Sweep** | Parameter sensitivity analysis. Sweeps any model parameter over a user-defined range, runs the full ODE at each point, and plots final titer, peak VCD, and metabolite endpoints. |
| **PC-dFBA** | Live SVG metabolic flux map (10 metabolite nodes, 16 reaction edges). Shows glycolysis → TCA → amino acid metabolism flux magnitudes and directions at any culture time point. PC trajectory visualisation included. |
| **GEM Red.** | Interactive GEM Reduction pipeline. Runs iCHOv1 through the full 5-step Antonakoudis & Richelle (2026) algorithm via FastAPI + COBRApy backend. Incremental step results with compression waterfall chart. |

---

## Scientific Background

### The ODE Kinetic Core (Eqs. 1–26)

Eight state variables evolve over time in a fed-batch bioreactor:

| Variable | Meaning |
|----------|---------|
| X_v | Viable cell density (10⁶ cells/mL) |
| X_d | Dead cell density |
| X' | Lysed cell density |
| Glc | Glucose concentration (mM) |
| Lac | Lactate concentration (mM) |
| Gln | Glutamine concentration (mM) |
| Glu | Glutamate concentration (mM) |
| NH₄⁺ | Ammonium concentration (mM) |

Glucose consumption is coupled to lactate production via a Warburg-like overflow stoichiometry that shifts as lactate accumulates. Glutamine spontaneous degradation to glutamate and ammonium (Tritsch & Moore 1986) is modelled explicitly.

### Net Growth Rate Models

The growth rate μ_net drives every metabolite flux. Three implementations are available:

1. **Sigmoid** — exact logistic basis function formulation from Eq. 14 of the paper
2. **Monod proxy** — biologically realistic Monod saturation + inhibition terms tuned to CHO DG44 fed-batch dynamics
3. **Surrogate Neural Network** — 5→8→1 MLP auto-calibrated via Adam-SGD at app startup, demonstrating data-driven hybrid modeling

### PC-dFBA (Eqs. 27–33)

Reduces the FBA flux space to principal components via PCA, then tracks the PC trajectory over the culture time course. Implemented analytically for the condensed 10-metabolite/16-reaction CHO network — no LP solver required in the browser.

### GEM Reduction (Antonakoudis & Richelle 2026)

| Step | Method | Description |
|------|--------|-------------|
| 0 | Model curation | Apply CHO-DG44 auxotrophies; set exchange bounds from 57 measured metabolites |
| 1 | Infeasibility resolution | Slack LP minimising total infeasibility across all time points |
| 2 | Exchange pruning | MILP: smallest subset of exchanges permitting feasible flux at every time point |
| 3 | Transport deduplication | Biological preference hierarchy eliminates redundant transport mechanisms |
| 4 | pFBA trimming | Remove zero-flux reactions across all parsimonious FBA solutions |
| 5 | Loop removal | Remove reactions active only under thermodynamically infeasible cycling |

**Result: 6,663 → 503 reactions (92.5% compression)**, retaining 105/155 metabolic tasks.

---

## Architecture

```
React 19 + Vite 7           ← ODE simulator, MetRaC, PC-dFBA, surrogate NN
     │                          (all run client-side, no server round-trips)
     │
     ├─ /gem/*  → Express proxy → FastAPI :8082
     │                              └─ COBRApy + HiGHS/GLPK
     │                                 (GEM reduction pipeline, ~22 s)
     └─ /api/*  → Express API   ← health, future: persistence, auth
```

**Client-side computation:** ODE integrator (RK4, dt=0.01 days), MetRaC Gaussian Process, PC-dFBA analytical solver, and surrogate NN all run in TypeScript in the browser. A 14-day simulation completes in <50 ms.

**Server-side COBRApy:** GEM reduction requires GLPK/HiGHS. Runs on FastAPI via a job-queue pattern: `POST /gem/run-pipeline` returns a `job_id`; the client polls `GET /gem/job/{job_id}` every 1.5 s; step results stream incrementally.

### Monorepo Structure

```
cho_digital_twin-/
├── artifacts/
│   ├── cho-digital-twin/     # React 19 + Vite frontend
│   │   └── src/
│   │       ├── lib/          # Core: simulator, ODE solver, models, FBA, MetRaC, NN
│   │       ├── pages/        # 7 application tabs
│   │       └── components/   # Radix UI + custom components
│   ├── api-server/           # Express.js proxy + API (TypeScript)
│   ├── gem-service/          # FastAPI + COBRApy GEM reduction (Python)
│   └── mockup-sandbox/       # Component development sandbox
├── lib/
│   ├── api-spec/             # OpenAPI schema + Orval codegen config
│   ├── api-zod/              # Generated Zod validation schemas
│   ├── api-client-react/     # React Query API client hooks
│   └── db/                   # Drizzle ORM schema (PostgreSQL)
├── scripts/                  # Utility scripts (tsx)
├── package.json              # pnpm workspace root
├── pnpm-workspace.yaml       # Workspace config + supply-chain security
└── tsconfig.base.json        # Shared TypeScript strict config
```

### Tech Stack

**Frontend**
- React 19.1 + Vite 7.3 + TypeScript 5.9
- Radix UI (45+ primitives) + Tailwind CSS 4 + Framer Motion
- Recharts 2.15 (simulation charts)
- TanStack React Query v5 + Zod v3
- React Hook Form + Wouter (routing)

**Backend**
- Express.js 5.2 (TypeScript proxy + API)
- FastAPI 0.136 + Uvicorn (Python GEM service)
- COBRApy 0.31 + HiGHS (MILP solver) + GLPK
- Pino (structured logging), http-proxy-middleware

**Build & Tooling**
- pnpm workspaces with supply-chain hardening (1-day minimum package release age)
- Orval (OpenAPI → TypeScript + Zod codegen)
- Drizzle ORM (PostgreSQL, future persistence layer)
- esbuild + Prettier

---

## Getting Started

### Prerequisites

- [Node.js](https://nodejs.org/) ≥ 20
- [pnpm](https://pnpm.io/) ≥ 9
- [Python](https://www.python.org/) ≥ 3.11

### Install dependencies

```bash
# JavaScript/TypeScript packages
pnpm install

# Python packages for the GEM service
cd artifacts/gem-service
pip install -r requirements.txt
cd ../..
```

### Run in development

```bash
# Start the React frontend + Express proxy
pnpm --filter ./artifacts/cho-digital-twin dev
pnpm --filter ./artifacts/api-server dev

# Start the FastAPI GEM service (required for the GEM Red. tab only)
cd artifacts/gem-service
uvicorn main:app --port 8082 --reload
```

The app will be available at `http://localhost:5173`.

> **Note:** The Simulator, MetRaC, PC-dFBA, Sweep, Equations, and Parameters tabs all run entirely client-side and work without the Python backend. Only the **GEM Red.** tab requires the FastAPI service.

### Build for production

```bash
pnpm build
```

---

## Scientific Citations

If you use this application in research, please cite the original papers:

> **Richelle A., Andersson D., Antonakoudis A., Jakobsson J., Pijeaud S., Vernersson A., Trygg J.** (2025). Integrating mechanistic and machine learning models for predictive digital twins of CHO cell culture. *bioRxiv* 2025.11.24.690194. https://doi.org/10.1101/2025.11.24.690194

> **Antonakoudis A., Richelle A.** (2026). Systematic data-driven genome-scale metabolic model reduction for bioprocess modeling: CHO culture case study. *npj Systems Biology and Applications*. https://doi.org/10.1038/s41540-026-00704-4

> **Hefzi H. et al.** (2016). A Consensus Genome-scale Reconstruction of Chinese Hamster Ovary Cell Metabolism. *Cell Systems* 3(5):434–443. https://doi.org/10.1016/j.cels.2016.10.020

---

## License

MIT License — Copyright (c) 2026 **Kemal Yaylali** ([support@yaylali.uk](mailto:support@yaylali.uk)) | ORCID: [0000-0003-1190-7807](https://orcid.org/0000-0003-1190-7807)

The application software is released under the MIT License. The underlying mathematical models, parameters, and algorithms are the intellectual work of the scientific paper authors cited above. See [LICENSE](LICENSE) for full terms.

---

## Author

**Kemal Yaylali**
- Email: [support@yaylali.uk](mailto:support@yaylali.uk)
- ORCID: [https://orcid.org/0000-0003-1190-7807](https://orcid.org/0000-0003-1190-7807)
- GitHub: [lynchaos](https://github.com/lynchaos)
