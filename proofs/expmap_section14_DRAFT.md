# DRAFT — not filed, awaiting Kit's approval (prepared 2026-10-15)

Proposed §14 for `proofs/sqrt_fw_clifford.md`, plus a companion
sketch for the 1820 (EPS) note §6. Sources: io-word sweep runs
815 (1820 side) and 816 (function side), ATTEMPTS.md entries of
2026-10-15 ~13:41–14:20 PDT; scripts run815_expmap_verify.py /
run816_expmap_partB.py and their results JSON in
`goals/r-theory-volume-i-audit/hidden_files/io_word_sweep/`.

---

## Addendum 2026-10-15 (VII) — the exponential map

### 14. exp(θJ) is sectorwise Euler (PROVED, operator-algebra level)

From the §2 premises — $J = S \cdot G$, $J^2 = -M_{FW}$,
$[J, M_{FW}] = 0$, $M_{FW}$ central — the exponential series
collapses by eigenspace of $M_{FW}$. Since
$J^{2k} = (-M_{FW})^k$ and $J^{2k+1} = J(-M_{FW})^k$:

$$\exp(\theta J) = A(\theta) + J\,B(\theta),$$

with $A, B$ diagonal: on the $FW = +1$ sectors
$(\cos\theta, \sin\theta)$ — rotations; on the $FW = -1$
sectors $(\cosh\theta, \sinh\theta)$ — boosts. One generator,
and the FlatWave sign selects whether its one-parameter group
is compact or noncompact: the §2 $\mathrm{Cl}(2,0)$ /
$\mathrm{Cl}(1,1)$ sector split at the group level.

**Block anatomy (exact, quad model).** $\pi$-transport pairs
quadrants of equal $FW$ sign, so $\exp(\theta J)$ splits: the
(I, III) pair block is $\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$
(det 1, unitary — a genuine rotation); the (II, IV) pair block
is $\begin{pmatrix} \cosh\theta & -i\sinh\theta \\ i\sinh\theta & \cosh\theta \end{pmatrix}$
(det 1, not unitary, and it exactly preserves the indefinite
Hermitian form $\eta = \mathrm{diag}(1, -1)$ in that pair
ordering, $E^\dagger \eta E = \eta$ — a Lorentz boost for
$\eta$, not merely cosh/sinh coefficients).

**Verification.** Premises re-verified exactly in both sympy
and Sage 9.5 (on toe): $S^2 = I$, $G^2 = M_{FW}$,
$\{S, G\} = 0$, $J^2 = -M_{FW}$, $[J, M_{FW}] = 0$.
Direct exact matrix exp (quad model) equals $A + JB$
entrywise; the ODE characterization $E(0) = I$, $E' = J \cdot E$
holds entrywise exactly in both the quad (4-quadrant) and oct
(8-octant) models, in sympy and Sage — the defining
characterization of $\exp(\theta J)$, so the oct model is exact,
not merely numeric. Numeric cross-checks: mpmath at 120 dps,
$\theta = 0.7$: max deviation exactly 0.0 in both models;
float64 scipy at $\theta \in \{-1.3, 0.7, 2.0\}$: max diff
$1.33 \times 10^{-15}$; group law
$\exp(0.3J)\exp(0.9J) = \exp(1.2J)$ at $1.11 \times 10^{-15}$.

Disclosure: sympy's plain `simplify` did not close 4 entries of
the direct-exp comparison to cos/cosh form (the same
presentation artifact as §§8/13); `rewrite(exp) + simplify`
closes every entry to exactly 0. The raw-simplify flag is
recorded False in the run's results file, not hidden. The
originating chat's $10^{-165}$ agreement figure was not
reproduced as such — the independent checks are the 120-dps
exact 0.0 and the float $10^{-15}$ figures above; no figure is
rounded up beyond what was measured.

Scope: PROVED at the operator-algebra level from the §2
premises (themselves exact-verified in-model); the finite
models are the filed instances of that algebra. No physical
identification is claimed. This is function/sector-model
algebra, not a physical input: the physical C8 identification,
$m_0/m_2$, and the o-ambiguity remain OPEN, untouched.

---

## Companion sketch — 1820 (EPS) note §6 (CHECKED)

The discrete root tower $\{I, \Gamma, M_\varepsilon, \Gamma^3\}$
on the 1820 is the integer-time sampling of a continuous flow.
With $K = \sum_k \kappa_k P_k$, $\kappa$ the principal angles
of $\gamma = (-\pi/2, 0, 0, \pi/2, \pi/2, \pi, \pi, -\pi/2)$:
spectrum $-\pi/2 \times 970$, $0 \times 600$, $\pi/2 \times 34$,
$\pi \times 216$ (from $\mathrm{tr}\,P_k = [554, 416, 184, 32, 2, 32, 184, 416]$);
$K$ Hermitian, $\mathrm{tr}\,K = -252\pi$. The flow
$U(t) = \exp(itK) = \sum_k e^{it\kappa_k} P_k$ satisfies
$U(1) = \Gamma$ ($1.1 \times 10^{-16}$),
$U(2) = M_\varepsilon$ ($2.5 \times 10^{-16}$),
$U(3) = \Gamma^3$ ($2.2 \times 10^{-15}$ via $M_\varepsilon \Gamma$),
$U(4) = I$ ($1.7 \times 10^{-15}$) — period exactly 4;
$U(1/2)^2 = \Gamma$ ($3.7 \times 10^{-15}$); group law and
unitarity at generic $t$ to $\sim 2.7 \times 10^{-15}$;
$\max|\mathrm{Im}\,U(2)| = 1.0$.

Scope: CHECKED on the archived float matrices (the $P_k$
premises are float-verified, so not PROVED). The only choice
in the construction is the log branch — per-sector $2\pi\mathbb Z$
windings, the §12 branch data as winding choices; principal
branch used. Physical C8 identification, $m_0/m_2$, and the
o-ambiguity remain OPEN, untouched.
