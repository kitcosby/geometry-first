# Exercise: Maxwell-Form Structures from the Core Operators

Date: 2026-10-12. Scope (per Kit's certified spine + algebra-first rule):
this is a REPRESENTATION WITNESS / formal exercise. Every operator
identity below is PROVED. The identification of operator families with
E and B is the exercise's declared dictionary (ASSERTED). Maxwell's
equations as spacetime PDEs are NOT derived — see the closing section.

## Step 1 (PROVED): the operators form a closed first-order system
In QI (t = tan(x/2), d/dx = ((1+t^2)/2) d/dt), the primitives obey
Riccati equations in two sign classes:
  srx' = -(1+srx^2)/2,   crx' = -(1+crx^2)/2   (cot type, decreasing)
  sxp' = +(1+sxp^2)/2,   cxp' = +(1+cxp^2)/2   (tan type, increasing)
and the UNA pair obeys log-derivative laws whose rates are exactly the
diagonal cross-family sums:
  urx'/urx = -(srx+crx)/2,   uxp'/uxp = +(cxp+sxp)/2.
(All verified symbolically.)

## Step 2 (PROVED): linearization = the wave (oscillator) equation
Riccati y' = -(1+y^2)/2 linearizes under y = 2 (log u)' to
  u'' + (1/4) u = 0,
the harmonic oscillator — the separated mode equation of the 1D wave
equation. Indeed srx = 2 (d/dx) log sin(x/2), a log-derivative of an
oscillator solution; the operators are Riccati variables of wave modes
at half-angle frequency.

## Step 3 (PROVED): curl eigenstructure / helicity
The so(3) placement cross matrix Y (omega = operator vector, as in the
uxp placement) has charpoly lambda (lambda^2 + theta^2),
theta^2 = |omega|^2: eigenvalues 0, +i theta, -i theta. The +/- i
eigenspaces are the two helicity (circular polarization) sectors of a
Fourier-mode curl operator. (Matches the certified witness
"Q^2+S^2=1 matches circular polarization.")

## Step 4 (PROVED identity): duality
The quarter-turn swaps the families: srx(x+pi/2) = crx(x),
cxp(x+pi/2) = sxp(x) (H4). A translation of the angle acts as a
discrete duality rotation exchanging the two families.

## Step 5 (ASSERTED dictionary): the Maxwell form
Source-free Maxwell for one Fourier mode (c=1): k x E = w B,
k x B = -w E, k.E = k.B = 0. Declare: family A <-> E, family B <-> B,
placement vector <-> k, quarter-turn <-> duality rotation. Under this
dictionary the operator system mirrors the Maxwell mode system:
wave modes (Step 2), helicity sectors (Step 3), duality (Step 4),
closed first-order dynamics (Step 1).

## What a derivation would still require (OPEN)
- A derived spacetime: the operators depend on ONE variable x; time
  and space (and the mode vector k) are imported by the ansatz.
- Field equations for all modes + sources (J) + gauge structure.
- A proved map from operator quantities to measured E, B.
None of these exist yet. Verdict: Maxwell's FORM is reproduced as a
witness under a declared dictionary; Maxwell's EQUATIONS are not
derived from the core operators. Identities PROVED; physics ASSERTED.
Script: /tmp/maxwell_ex.py (this run).
