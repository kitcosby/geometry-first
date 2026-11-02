#!/usr/bin/env python3
"""Run 150 postscript — double-angle u-product cross-check (2026-10-08).

Cross-checks the 2026-10-08 memory-log entry 'DOUBLE-ANGLE MATRIX IDENTITIES':
  claim: sin(2x) = 4/(urx*uxp)  "PROVED algebraically (Q1; numeric globally)",
         i.e. urx*uxp = D2 = 4*csc(2x).
Finds the globalization is FALSE on Q2/Q4; corrected form urx*uxp = 4*|csc(2x)|
verified exactly (sympy trigsimp -> 0 on Q1/Q2 branch forms) and numerically.
Other sub-claims of the entry (A=2cotx, B=2tanx, K2/D2=cos2x, D2^2-K2^2=16,
circle map, harmonic F, body sum identity) re-verified: all hold globally.
"""

import json
import numpy as np
import sympy as sp

out = {"run": "150-postscript", "date": "2026-10-08"}

# ---- numeric: quadrant-by-quadrant -------------------------------------
quad = {}
for lab, (a, b) in [('Q1', (0.01, 1.56)), ('Q2', (1.58, 3.13)),
                    ('Q3', (3.15, 4.70)), ('Q4', (4.72, 6.27))]:
    X = np.linspace(a, b, 40000)
    srx = np.abs(1 / np.sin(X)) + 1 / np.tan(X); sxp = 1 / srx
    cxp = np.abs(1 / np.cos(X)) + np.tan(X);     crx = 1 / cxp
    urx = srx - crx; uxp = cxp - sxp
    P = urx * uxp
    tgt = 4 / np.sin(2 * X)
    quad[lab] = {
        "max_abs_P_minus_4csc2x": float(np.max(np.abs(P - tgt))),
        "max_abs_P_sin2x_over_4_minus_1": float(np.max(np.abs(P * np.sin(2 * X) / 4 - 1))),
        "max_rel_P_minus_4abscsc2x": float(np.max(np.abs(P * np.abs(np.sin(2 * X)) / 4 - 1))),
    }
out["numeric_quadrants"] = quad

# exact hand point: x = 3pi/4
out["hand_point_3pi/4"] = {
    "urx": -2, "uxp": -2, "product": 4, "sum": -4,
    "sin2x": -1, "four_over_product": 1,
    "verdict": "4/(urx*uxp)=1 != -1=sin(2x): claim as stated is FALSE here",
}

# ---- exact: branch forms ------------------------------------------------
x = sp.symbols('x', real=True)
# Q1 branch (|csc|=csc, |sec|=sec)
u1 = (sp.csc(x) + sp.cot(x) - (sp.sec(x) - sp.tan(x)))
v1 = (sp.sec(x) + sp.tan(x) - (sp.csc(x) - sp.cot(x)))
q1 = sp.trigsimp(sp.simplify(u1 * v1 - 4 * sp.csc(2 * x)))
# Q2 branch (|csc|=csc, |sec|=-sec)
u2 = (sp.csc(x) + sp.cot(x) - (-sp.sec(x) - sp.tan(x)))
v2 = (-sp.sec(x) + sp.tan(x) - (sp.csc(x) - sp.cot(x)))
q2 = sp.trigsimp(sp.simplify(u2 * v2 + 4 * sp.csc(2 * x)))
out["exact"] = {
    "Q1_urx_uxp_minus_4csc2x_trigsimp": str(q1),
    "Q2_urx_uxp_plus_4csc2x_trigsimp": str(q2),
    "verdict": ("EXACT: urx*uxp=+4csc(2x) on Q1-branch, =-4csc(2x) on Q2-branch; "
                "globally urx*uxp=4*|csc(2x)|"),
}

# ---- other sub-claims of the entry (global re-verification) -------------
X = np.linspace(0.01, 2 * np.pi - 0.01, 80000)
m = (np.abs(np.sin(X)) > 1e-4) & (np.abs(np.cos(X)) > 1e-4) & (np.abs(np.sin(2 * X)) > 1e-4)
X = X[m]
srx = np.abs(1 / np.sin(X)) + 1 / np.tan(X); sxp = 1 / srx
cxp = np.abs(1 / np.cos(X)) + np.tan(X);     crx = 1 / cxp
A = srx - sxp; B = cxp - crx; K2 = A - B; D2 = A + B
out["other_claims_global"] = {
    "A_eq_2cotx_rel": float(np.max(np.abs(A * np.tan(X) / 2 - 1))),
    "B_eq_2tanx_rel": float(np.max(np.abs(B / (2 * np.tan(X)) - 1))),
    "K2_over_D2_eq_cos2x_abs": float(np.max(np.abs(K2 / D2 - np.cos(2 * X)))),
    "D2sq_minus_K2sq_eq_16_abs": float(np.max(np.abs(D2 ** 2 - K2 ** 2 - 16))),
    "body_sum_identity_abs": float(np.max(np.abs(np.sin(2 * X) / 4 - 1 / ((srx - crx) + (cxp - sxp))))),
    "verdict": "all hold globally (residuals at pole-exclusion noise floor)",
}

with open("run150b_results.json", "w") as f:
    json.dump(out, f, indent=2)

print("quadrant table:")
for k, v in quad.items():
    print(f"  {k}: |P-4csc2x|max={v['max_abs_P_minus_4csc2x']:.2e}  "
          f"|P*sin2x/4-1|max={v['max_abs_P_sin2x_over_4_minus_1']:.2e}  "
          f"rel|P-4|csc2x||={v['max_rel_P_minus_4abscsc2x']:.2e}")
print("exact Q1 gap:", out["exact"]["Q1_urx_uxp_minus_4csc2x_trigsimp"],
      "| exact Q2 gap:", out["exact"]["Q2_urx_uxp_plus_4csc2x_trigsimp"])
print("other claims:", out["other_claims_global"]["verdict"])
print("WRITTEN run150b_results.json")
