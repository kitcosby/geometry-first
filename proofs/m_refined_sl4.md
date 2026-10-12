# M_refined(x): gl(4) raw, sl(4,R) traceless-corrected

Date: 2026-10-12. Status: PROVED (Sage 9.5 exact: QQ(t) for det,
QQbar for closure/Killing). Matrix supplied by Kit (image upload).

M_refined(x) =
[ sin x,   cxp,    -sxp,    uxp   ]
[ srx,     cos x,   cos 2x, sawup ]   (1/urx = sawup, 1/uxp = sawdown)
[ -crx,    sin 2x,  0,      0     ]
[ urx,     sawdown, 0,      0     ]

## Results
- tr M = sin x + cos x (not traceless).
- det M(t) = -15 (t + 5/3)(t^2 - 1/5) / [ (t+1)(t^2+1)^2 ], t = tan(x/2),
  an exact rational identity over QQ(t). It factors over the zero block:
  det = A * B with A = -sxp*sawup - uxp*cos(2x) (upper-right block det),
  numerator -5(t^2 - 1/5)/[(t-1)(t^2+1)]; B = -crx*sawdown - sin(2x)*urx
  (lower-left block det), numerator 3(t-1)(t+5/3)/[(t+1)(t^2+1)].
  In QI, det = 0 at t^2 = 1/5, i.e. x = 2 arctan(1/sqrt5) ~ 48.19 deg
  (the upper block degenerates). B also vanishes at t = -5/3 (negative-t
  quadrants).
- Samples (5 exact QI points over QQbar): commutator closure dim 16,
  closed => the raw matrix generates gl(4,R).
- Traceless correction M0 = M - ((sin x + cos x)/4) I: closure dim 15,
  closed, Killing signature (9, 6, 0) => sl(4,R) PROVED (among A3 real
  forms: sl(4,R) (9,6), su(2,2) (8,7), su(3,1) (6,9), su(4) (0,15)).

## Disclosures
- First attempt in sympy produced no output in 20+ min (nested-radical
  rank calls) and was stopped; all results here are Sage.
- Sage over AA failed at span() with the known toe numpy binary
  incompatibility; QQbar path used instead. Polynomial .roots() is
  unsupported over AA; signature read from QQbar eigenvalues with exact
  real-part signs, multiplicities included.
Scripts: /tmp/mref_sage.py (det), /tmp/mref_sage2.py (closures),
/tmp/mref_sage4.py (Killing), this run.
