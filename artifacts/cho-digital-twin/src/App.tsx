
import { useState } from "react";
import SimulatorPage from "@/pages/SimulatorPage";
import EquationsPage from "@/pages/EquationsPage";
import ParametersPage from "@/pages/ParametersPage";
import MetRaCPage from "@/pages/MetRaCPage";
import SweepPage from "@/pages/SweepPage";
import PcdFBAPage from "@/pages/PcdFBAPage";

type Tab = "simulator" | "equations" | "parameters" | "metrac" | "sweep" | "pcdfba" | "about";

const TABS: { id: Tab; label: string; icon: string }[] = [
  { id: "simulator",  label: "Simulator",  icon: "⚙" },
  { id: "equations",  label: "Equations",  icon: "∑" },
  { id: "parameters", label: "Parameters", icon: "⊞" },
  { id: "metrac",     label: "MetRaC",     icon: "≈" },
  { id: "sweep",      label: "Sweep",      icon: "⊹" },
  { id: "pcdfba",     label: "PC-dFBA",    icon: "⬡" },
  { id: "about",      label: "About",      icon: "ℹ" },
];

function AboutPage() {
  return (
    <div className="about-page">
      <h1>CHO Cell Culture Digital Twin</h1>
      <p className="about-tagline">
        Replication of the hybrid modeling framework for predictive digital twins of CHO cell culture
      </p>

      <div className="about-ref">
        <strong>Source paper:</strong>{" "}
        Richelle A., Andersson D., Antonakoudis A., Jakobsson J., Pijeaud S., Vernersson A., Trygg J. (2025)
        <em> A Hybrid Modeling Framework for Predictive Digital Twins of CHO Cell Culture.</em>{" "}
        bioRxiv 2025.11.24.690194. doi:{" "}
        <a href="https://doi.org/10.1101/2025.11.24.690194" target="_blank" rel="noopener">
          10.1101/2025.11.24.690194
        </a>
      </div>

      <h2>What is replicated</h2>
      <table className="about-table">
        <thead><tr><th>Component</th><th>Equations</th><th>Status</th></tr></thead>
        <tbody>
          <tr><td>ODE biomass population model</td><td>Eqs. 8–14</td>
            <td className="status-yes">✓ Exact Table 1 parameters; B_max cap added</td></tr>
          <tr><td>ODE FLEX metabolite model</td><td>Eqs. 15–26</td>
            <td className="status-yes">✓ Exact Table 1 parameters</td></tr>
          <tr><td>Gln chemical degradation</td><td>Eqs. 17, 18</td>
            <td className="status-yes">✓ k_Gln_deg = 0.006 day⁻¹ (Table 1); NH₄⁺ co-production</td></tr>
          <tr><td>Gln → Glu transamination coupling</td><td>Eqs. 17–18, 24–25</td>
            <td className="status-yes">✓ Y_Glu_Gln = 0.67 (Table 1); Glu rise on Gln depletion</td></tr>
          <tr><td>RK4 integrator</td><td>—</td>
            <td className="status-yes">✓ dt = 0.005 days</td></tr>
          <tr><td>Fed-batch mass balance with bolus feeds</td><td>Eq. 34</td>
            <td className="status-yes">✓ Editable bolus schedule; exact volume mixing</td></tr>
          <tr><td>Luedeking-Piret product model</td><td>§2.4 extension</td>
            <td className="status-yes">✓ dTit/dt = (α·μ_net + β)·Xv — growth + non-growth associated</td></tr>
          <tr><td>μ_net growth rate (§2.3)</td><td>§2.3</td>
            <td className="status-partial">⚠ Three modes: Sigmoid · Monod proxy · Surrogate NN (auto-calibrated from Monod proxy; full NN needs proprietary 23-batch dataset)</td></tr>
          <tr><td>MetRaC rate estimation + q_p</td><td>§2.2</td>
            <td className="status-partial">⚠ Two methods: (a) kernel-smooth finite-diff (fast); (b) SE-kernel GP with analytical derivative posterior — proper Bayesian CIs, length-scale via marginal likelihood. Full nested-sampling B-splines not implemented.</td></tr>
          <tr><td>PC-dFBA hybrid LP</td><td>Eqs. 27–33</td>
            <td className="status-partial">⚠ Condensed 10-met/16-rxn network; analytical mass-balance solver; PC analysis from flux trajectories (iCHO2441 + full NN loadings not available)</td></tr>
        </tbody>
      </table>

      <h2>Framework overview</h2>
      <div className="framework-diagram">
        <div className="fw-node fw-input">
          Experimental data<br/><small>23 CHO fed-batch runs</small>
        </div>
        <div className="fw-arrow">→</div>
        <div className="fw-node fw-proc">
          MetRaC<br/><small>§2.2 Bayesian rates</small>
        </div>
        <div className="fw-arrow">→</div>
        <div className="fw-col">
          <div className="fw-node fw-model">ODE Biomass<br/><small>Eqs. 8–14 ✓</small></div>
          <div className="fw-node fw-model">ODE FLEX<br/><small>Eqs. 15–26 ✓</small></div>
          <div className="fw-node fw-model">VCD NN (μ_net)<br/><small>§2.3 → Monod proxy ⚠</small></div>
          <div className="fw-node fw-model">PC-dFBA<br/><small>Eqs. 27–33 ✗</small></div>
        </div>
        <div className="fw-arrow">→</div>
        <div className="fw-node fw-output">
          Digital Twin<br/><small>VCD · titer · metabolites</small>
        </div>
      </div>

      <h2>Implementation notes</h2>
      <ul className="about-list">
        <li>
          <strong>Gln degradation &amp; transamination:</strong> Glutamine degrades abiotically at k_Gln_deg = 0.006 day⁻¹
          (Table 1, Eqs. 17–18), producing NH₄⁺. The Y_Glu_Gln = 0.67 stoichiometric coupling drives Glu accumulation
          as Gln depletes, clearly visible in the Gln/Glu chart after day 8.
        </li>
        <li>
          <strong>B_max cap:</strong> Biomaterial B is capped at B_max (default 500) inside deathRate(),
          preventing the kd1·B term from growing unboundedly and producing unrealistic death rates after day 10.
        </li>
        <li>
          <strong>Luedeking-Piret product model:</strong> dTit/dt = (α·μ_net + β)·Xv, distinguishing
          growth-associated (α = q_p_growth) from non-growth-associated (β = q_p) productivity.
          The product-specific rate q_p is also tracked as a state in simulator output for MetRaC comparison.
        </li>
        <li>
          <strong>Three μ_net modes:</strong> (1) <em>Sigmoid</em> — bare Gaussian-sum baseline, no nutrient feedback;
          (2) <em>Monod proxy</em> — μ_sigmoid × Monod(Glc) × Monod(Gln) × Inhibition(Lac) × Inhibition(NH₄⁺),
          gives biologically realistic VCD without training data;
          (3) <em>Surrogate NN</em> — MLP auto-calibrated from the Monod proxy teacher signal, demonstrating
          the §2.3 architecture.
        </li>
        <li>
          <strong>Feed schedule editor:</strong> The Simulator sidebar includes an editable bolus table —
          add, remove, or reset glucose/glutamine feeds and observe the effect on metabolite trajectories in real time.
        </li>
        <li>
          <strong>Reference run:</strong> Click <em>📌 Set Reference</em> to pin the current run as a ghost
          overlay on all charts; change parameters and compare side-by-side.
        </li>
        <li>
          <strong>Parameter sweep:</strong> The Sweep tab scans any single Table 1 parameter across a range
          and plots a chosen output metric (peak VCD, final titer, etc.) against it — useful for sensitivity analysis.
        </li>
        <li>
          <strong>MetRaC q_p estimation:</strong> Enable Titer CV &gt; 0 in the MetRaC noise panel to reveal
          the q_p rate chart — MetRaC estimates the product-specific rate from noisy titer measurements
          and overlays the Luedeking-Piret ODE truth.
        </li>
        <li>
          <strong>JSON export / import:</strong> Use the Export/Import buttons at the bottom of the Simulator
          sidebar to serialise the full parameter and bolus state to clipboard (JSON), then restore it later.
        </li>
        <li>
          <strong>CSV export:</strong> Click ↓ CSV in the Simulator to download all state variables
          and specific rates at every output step.
        </li>
        <li><strong>Overflow metabolism:</strong> Aerobic lactate production when glucose uptake exceeds oxidative capacity (Eqs. 21–23)</li>
        <li><strong>Lactate switch:</strong> Cells re-consume lactate when glucose is limiting (Eq. 22)</li>
        <li><strong>Biomaterial inhibition:</strong> Accumulating by-products increase death and lysis (Eqs. 12–13)</li>
      </ul>

      <h2>What still needs experimental data</h2>
      <ul className="about-list">
        <li>
          <strong>§2.3 NN weights:</strong> The neural network that predicts μ_net from specific rates and metabolite
          concentrations was trained on 23 fed-batch CHO runs (AstraZeneca/Sartorius proprietary dataset).
          The Monod nutrient coupling is a structural substitute.
        </li>
        <li>
          <strong>PC-dFBA (Eqs. 27–33):</strong> Requires the CHO genome-scale model (e.g., iCHO2441, ~2000 reactions),
          a linear programming solver, PCA loadings predicted by a second NN, and MetRaC-derived exchange rates
          as boundary conditions.
        </li>
      </ul>
    </div>
  );
}

export default function App() {
  const [activeTab, setActiveTab] = useState<Tab>("simulator");

  return (
    <div className="app-shell">
      <header className="app-header">
        <div className="app-header-inner">
          <div className="app-logo">
            <span className="app-logo-icon">🔬</span>
            <span className="app-logo-text">
              CHO Digital Twin
              <span className="app-logo-sub"> — Richelle et al. (2025) Replication</span>
            </span>
          </div>
          <nav className="app-nav">
            {TABS.map((tab) => (
              <button key={tab.id}
                className={`nav-tab ${activeTab === tab.id ? "nav-tab-active" : ""}`}
                onClick={() => setActiveTab(tab.id)}>
                <span className="nav-tab-icon">{tab.icon}</span>
                {tab.label}
              </button>
            ))}
          </nav>
        </div>
      </header>

      <div className="app-content">
        {activeTab === "simulator"  && <SimulatorPage />}
        {activeTab === "equations"  && <EquationsPage />}
        {activeTab === "parameters" && <ParametersPage />}
        {activeTab === "metrac"     && <MetRaCPage />}
        {activeTab === "sweep"      && <SweepPage />}
        {activeTab === "pcdfba"     && <PcdFBAPage />}
        {activeTab === "about"      && <AboutPage />}
      </div>
    </div>
  );
}
