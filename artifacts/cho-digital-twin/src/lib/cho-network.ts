/**
 * Condensed CHO central carbon metabolism network
 * ~10 intracellular metabolites, 16 reactions (4 exchange + 12 internal)
 *
 * Metabolite/reaction set covers: glycolysis (lumped), pyruvate metabolism,
 * TCA cycle, malic enzyme anaplerosis, pyruvate carboxylase, and
 * glutamine/glutamate catabolism (GLS + GLUD).
 *
 * Sign conventions for exchange fluxes (match TimePoint in simulator.ts):
 *   EX_GLC > 0  = glucose uptake rate  (matches q_Glc positive = consumed)
 *   EX_LAC > 0  = net lactate secretion (matches q_Lac positive = produced)
 *   EX_GLN > 0  = glutamine uptake rate (matches q_Gln positive = consumed)
 *   EX_GLU > 0  = glutamate UPTAKE rate (matches q_Glu positive = consumed from broth)
 * Note: GLU produced via transamination is not captured in this condensed model.
 */

// ── Metabolite IDs ────────────────────────────────────────────────────────────
export const MET_IDS = [
  "PYR","LAC","AcCoA","OAA","CIT","AKG","SUCC","MAL","GLU","GLN",
] as const;
export type MetId = (typeof MET_IDS)[number];

export interface MetaboliteNode {
  id: MetId | string;
  label: string;
  x: number; y: number;          // SVG canvas position (centred on node)
  extracellular?: boolean;
}

// ── Reaction IDs ─────────────────────────────────────────────────────────────
export const RXNS_EXCHANGE = ["EX_GLC","EX_LAC","EX_GLN","EX_GLU"] as const;
export const RXNS_INTERNAL = [
  "LDH_f","LDH_r","PDH","CS","ICDH","AKGDH","SDH_FH","MDH","ME","PCX","GLS","GLUD",
] as const;
export const RXN_IDS = [...RXNS_EXCHANGE, ...RXNS_INTERNAL] as const;
export type RxnId = (typeof RXN_IDS)[number];

export const N_METS = MET_IDS.length;    // 10
export const N_RXNS = RXN_IDS.length;    // 16
export const N_EXC  = RXNS_EXCHANGE.length; // 4
export const N_INT  = RXNS_INTERNAL.length; // 12

// ── Reaction metadata ─────────────────────────────────────────────────────────
export interface ReactionMeta {
  id: RxnId;
  label: string;
  formula: string;
  isExchange: boolean;
  lb: number; ub: number;        // default flux bounds
  fromId: string; toId: string;  // for SVG arrow
}

export const REACTION_META: Record<RxnId, ReactionMeta> = {
  EX_GLC:  { id:"EX_GLC",  label:"EX_GLC",   formula:"Glc → 2 Pyr (glycolysis)",          isExchange:true,  lb:0,     ub:100, fromId:"GLC_EXT", toId:"PYR"     },
  EX_LAC:  { id:"EX_LAC",  label:"EX_LAC",   formula:"Lac ↔ (extracellular)",              isExchange:true,  lb:-100,  ub:100, fromId:"LAC",     toId:"LAC_EXT" },
  EX_GLN:  { id:"EX_GLN",  label:"EX_GLN",   formula:"Gln → (uptake)",                    isExchange:true,  lb:0,     ub:100, fromId:"GLN_EXT", toId:"GLN"     },
  EX_GLU:  { id:"EX_GLU",  label:"EX_GLU",   formula:"Glu (uptake from broth)",           isExchange:true,  lb:0,     ub:100, fromId:"GLU_EXT", toId:"GLU"     },
  LDH_f:   { id:"LDH_f",   label:"LDH →",    formula:"Pyr → Lac",                          isExchange:false, lb:0,     ub:100, fromId:"PYR",     toId:"LAC"     },
  LDH_r:   { id:"LDH_r",   label:"LDH ←",    formula:"Lac → Pyr",                          isExchange:false, lb:0,     ub:100, fromId:"LAC",     toId:"PYR"     },
  PDH:     { id:"PDH",     label:"PDH",       formula:"Pyr → AcCoA + CO₂",                 isExchange:false, lb:0,     ub:100, fromId:"PYR",     toId:"AcCoA"   },
  CS:      { id:"CS",      label:"CS",        formula:"OAA + AcCoA → Cit",                 isExchange:false, lb:0,     ub:100, fromId:"AcCoA",   toId:"CIT"     },
  ICDH:    { id:"ICDH",    label:"ICDH",      formula:"Cit → αKG + CO₂",                  isExchange:false, lb:0,     ub:100, fromId:"CIT",     toId:"AKG"     },
  AKGDH:   { id:"AKGDH",   label:"αKGDH",     formula:"αKG → Succ + CO₂",                 isExchange:false, lb:0,     ub:100, fromId:"AKG",     toId:"SUCC"    },
  SDH_FH:  { id:"SDH_FH",  label:"SDH/FH",    formula:"Succ → Mal",                        isExchange:false, lb:0,     ub:100, fromId:"SUCC",    toId:"MAL"     },
  MDH:     { id:"MDH",     label:"MDH",       formula:"Mal → OAA",                         isExchange:false, lb:0,     ub:100, fromId:"MAL",     toId:"OAA"     },
  ME:      { id:"ME",      label:"ME",        formula:"Mal → Pyr + CO₂",                   isExchange:false, lb:0,     ub:100, fromId:"MAL",     toId:"PYR"     },
  PCX:     { id:"PCX",     label:"PCX",       formula:"Pyr + CO₂ → OAA",                  isExchange:false, lb:0,     ub:100, fromId:"PYR",     toId:"OAA"     },
  GLS:     { id:"GLS",     label:"GLS",       formula:"Gln → Glu + NH₃",                   isExchange:false, lb:0,     ub:100, fromId:"GLN",     toId:"GLU"     },
  GLUD:    { id:"GLUD",    label:"GLUD",      formula:"Glu → αKG + NH₃",                  isExchange:false, lb:0,     ub:100, fromId:"GLU",     toId:"AKG"     },
};

// ── Stoichiometric matrix S (MET_IDS × RXN_IDS) ──────────────────────────────
// S[i][j] = stoichiometric coefficient of metabolite i in reaction j
// Positive = produced, Negative = consumed
//
//          EXgc EXlc EXgn EXgu  LDHf LDHr PDH  CS  ICDH AKGDH SDH  MDH  ME  PCX  GLS  GLUD
//          [0]  [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]  [10] [11] [12] [13] [14] [15]
export const S_MATRIX: number[][] = [
  //PYR:
  [ +2,  0,  0,  0,  -1, +1, -1,  0,  0,  0,  0,  0, +1, -1,  0,  0 ],
  //LAC:
  [  0, -1,  0,  0,  +1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  //AcCoA:
  [  0,  0,  0,  0,   0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0 ],
  //OAA:
  [  0,  0,  0,  0,   0,  0,  0, -1,  0,  0,  0, +1,  0, +1,  0,  0 ],
  //CIT:
  [  0,  0,  0,  0,   0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0 ],
  //AKG:
  [  0,  0,  0,  0,   0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0, +1 ],
  //SUCC:
  [  0,  0,  0,  0,   0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0 ],
  //MAL:   SDH_FH produces MAL (+1); AKGDH produces SUCC (0 here)
  [  0,  0,  0,  0,   0,  0,  0,  0,  0,  0, +1, -1, -1,  0,  0,  0 ],
  //GLU:   EX_GLU is +1 (uptake adds to intracellular pool)
  [  0,  0,  0, +1,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1 ],
  //GLN:
  [  0,  0, +1,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0 ],
];

// ── SVG Layout for metabolic map ──────────────────────────────────────────────
// Canvas: 760 × 400 px
export const SVG_METABOLITES: MetaboliteNode[] = [
  // Glycolysis / pyruvate hub
  { id:"GLC_EXT", label:"Glucose",   x: 70,  y: 60,  extracellular:true  },
  { id:"PYR",     label:"Pyruvate",  x: 220, y: 160                       },
  { id:"LAC",     label:"Lactate",   x: 100, y: 280                       },
  { id:"LAC_EXT", label:"Lactate",   x:  50, y: 375, extracellular:true  },
  { id:"AcCoA",   label:"Acetyl-CoA",x: 330, y: 280                       },

  // TCA cycle (clockwise: OAA→CIT→AKG→SUCC→MAL→OAA)
  { id:"OAA",     label:"OAA",       x: 430, y: 140                       },
  { id:"CIT",     label:"Citrate",   x: 570, y: 140                       },
  { id:"AKG",     label:"α-KG",      x: 590, y: 280                       },
  { id:"SUCC",    label:"Succinate", x: 500, y: 380                       },
  { id:"MAL",     label:"Malate",    x: 370, y: 380                       },

  // Amino acid arm
  { id:"GLN_EXT", label:"Glutamine", x: 760, y: 60,  extracellular:true  },
  { id:"GLN",     label:"Glutamine", x: 720, y: 160                       },
  { id:"GLU",     label:"Glutamate", x: 720, y: 280                       },
  { id:"GLU_EXT", label:"Glutamate", x: 790, y: 375, extracellular:true  },
];
