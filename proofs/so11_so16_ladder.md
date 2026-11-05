# The orthogonal ladder so(11)–so(16): path algebras and spinors

Extension of the so campaign (so(3)–so(10), see `so10_path.md`,
`spin10_spinor.md`, `so6_so5.md`, `so7_so4_so3.md`, `so9_spinor.md`) to
n = 11..16, by Kit's order of 2026-10-16 ("extend the SO ladder to follow
the generators all the way to 16"). All computations ran in Sage 9.5 on the
laptop (WSL `lampy`) — the campaign's engine — except where a section says
otherwise (the spinor work uses an exact monomial engine, cross-validated
against dense Sage; see §3). LOCAL FILE — not uploaded.

Constructions are Kit's, identical to the earlier rungs:

- **Path algebra.** $Y_i = f_i\,(E_{i,i+1} - E_{i+1,i})$, $i = 1..n-1$, with
  the 9 established weights at $x = 30^\circ$ cycled,
  $f_i = W[(i-1) \bmod 9]$, $W = (srx, sxp, cxp, crx, urx, sawup, uxp,
  sawdown, FW)$. At $30^\circ$: $srx = 2+\sqrt3$, $sxp = 2-\sqrt3$,
  $cxp = \sqrt3$, $crx = \sqrt3/3$, $urx = 2 + 2\sqrt3/3$,
  $uxp = 2\sqrt3 - 2$, $sawup = 1/urx$, $sawdown = 1/uxp$, $FW = 1$.
  (Kit's length-generic induction already proves the path generates
  $\mathfrak{so}(n)$ for every $n$; the Sage closures below are the
  campaign-style exact cross-check at the cycled weights.)
- **Spinor.** Gamma matrices by the campaign's Pauli recursion
  ($\gamma_{2j-1} = I^{\otimes(j-1)} \otimes \sigma_1 \otimes
  \sigma_3^{\otimes(m-j)}$, $\gamma_{2j} = I^{\otimes(j-1)} \otimes \sigma_2
  \otimes \sigma_3^{\otimes(m-j)}$ on $m$ qubits),
  $\Sigma_{ab} = \frac14[\gamma_a, \gamma_b] = \frac12 \gamma_a\gamma_b$
  ($a<b$), module dimension $2^m$, $m = \lfloor (n+1)/2 \rfloor$ in the
  campaign's convention (odd $n$: first $n$ gammas of $Cl(n+1)$, the volume
  element splits the module into two equal blocks, each carrying the full
  $\mathfrak{so}(n)$).

## 1. Path algebras — PROVED (Sage-exact over QQbar)

Lie closure by incremental RREF span over the antisymmetric coordinates,
exactly as in `so10_path.md`. The $n-1$ path generators are linearly
independent at every rung (gen-span $= n-1$), and the closure reaches the
full $\mathfrak{so}(n)$ dimension at every rung:

| n | gen-span | closure | target n(n−1)/2 | match | wall time |
|---|---|---|---|---|---|
| 11 | 10 | **55** | 55 | ✓ | 111.9 s |
| 12 | 11 | **66** | 66 | ✓ | 211.5 s |
| 13 | 12 | **78** | 78 | ✓ | 386.6 s |
| 14 | 13 | **91** | 91 | ✓ | 673.4 s |
| 15 | 14 | **105** | 105 | ✓ | 1360.8 s |
| 16 | 15 | **120** | 120 | ✓ | 1812.2 s |

Script: `so11_so16_path.py`.

## 2. Spinors — per-rung results

Checks run at every rung (labels in §3): gamma anticommutation exact;
$\Sigma_{ab}^2 = -\frac14 I$ for every pair $a<b$ and $\Sigma_{ab}$
anti-Hermitian; the Lie homomorphism
$[\Sigma_{ab}, \Sigma_{cd}] = \Sigma_{[e_{ab}, e_{cd}]}$ on **all**
unordered distinct basis pairs, zero failures; $\operatorname{span}
\{\Sigma_{ab}\}$ of dimension $\dim \mathfrak{so}(n)$; chirality/volume
element $\Gamma$ with $\Gamma^2 = I$, $\operatorname{tr} \Gamma = 0$,
$[\Gamma, \Sigma_{ab}] = 0$ for all pairs, and eigenspace split as listed.

| n | module | Σ rank | hom pairs (fails) | Γ split | Casimir vector | Casimir spinor |
|---|---|---|---|---|---|---|
| 11 | 64 | 55 | 1485 (0) | 32+32 | −10 | −55/4 |
| 12 | 64 | 66 | 2145 (0) | 32+32 | −11 | −33/2 |
| 13 | 128 | 78 | 3003 (0) | 64+64 | −12 | −39/2 |
| 14 | 128 | 91 | 4095 (0) | 64+64 | −13 | −91/4 |
| 15 | 256 | 105 | 5460 (0) | 128+128 | −14 | −105/4 |
| 16 | 256 | 120 | 7140 (0) | 128+128 | −15 | −30 |

Casimir normalization is the campaign's: basis $e_{ab} = E_{ab} - E_{ba}$,
$C_2 = \sum_{a<b} \rho(e_{ab})^2$. Vector values $-(n-1)\,I$ were computed
**exactly over QQ in Sage** at every rung. Spinor values
$-n(n-1)/8$ follow exactly from $\Sigma_{ab}^2 = -\frac14 I$ summed over
the $n(n-1)/2$ basis elements. Vector and spinor Casimirs agree only at
$n = 8$ (both $-7$), as in the earlier rungs; at $n = 16$ they are $-15$
vs $-30$.

### Invariant bilinear forms

On the full module an intertwiner $J$ (a product of gamma matrices over a
subset $S$ of indices, with the volume-element adjustment at odd $n$) was
found at every rung and verified exactly:
$J\,\Sigma_{ab} + \Sigma_{ab}^{T} J = 0$ for **all** pairs $a<b$.
$\tau$ is the commutation-pattern sign returned by the engine's
`tau_test`: the single sign with $J\gamma_c = s_c\,\gamma_c J$ and
$s_c = \tau \cdot \varepsilon_c$ for all $c$, where
$\varepsilon_c = +1$ for odd $c$, $-1$ for even $c$; it is NOT the
transpose sign. The transpose sign
$J^T = tJ$ is carried by the type column: $t = -1$ for the
antisymmetric (symplectic) rows and $t = +1$ for the symmetric rows
— so $t = (-1, -1, -1, +1, +1, +1)$ for $n = (11, \dots, 16)$,
differing from $\tau$ at $n = 12, 13, 14, 15$. (Correction
2026-10-16, verified by re-running `so11_so16_spinor.py`: the table
previously defined its $\sigma$ column as the transpose sign while
listing the engine's $\tau$ values.)

| n | S (indices of γ factors) | τ | type | block behaviour |
|---|---|---|---|---|
| 11 | {2,4,6,8,10} | −1 | antisymmetric (symplectic) | preserves blocks |
| 12 | {2,4,…,12} | +1 | antisymmetric (symplectic) | preserves blocks |
| 13 | {2,4,…,12} | +1 | antisymmetric (symplectic) | preserves blocks |
| 14 | {2,4,…,14} | −1 | symmetric | **swaps blocks** |
| 15 | {2,4,…,14} | −1 | symmetric | preserves blocks |
| 16 | {2,4,…,16} | +1 | symmetric | preserves blocks |

(Control: the same search at $n = 10$ returns $S = \{2,\dots,10\}$,
$\tau = -1$, antisymmetric, swaps blocks — the campaign's so(10) form;
at $n = 10$ the $\tau$ and transpose signs coincide, both $-1$.)

Block-level classification (the campaign's per-block convention; a
preserving $J$ is invertible, so its restriction to each block is a
nonzero invariant form of the same type):

- **n = 11** (block 32): symplectic; **unique** — solution-space solve
  (§3): symmetric dimension 0, antisymmetric dimension 1, at two primes.
- **n = 12** (blocks 32+32): symplectic on **each** block; **unique on
  each** — both blocks solved: symmetric 0, antisymmetric 1, two primes.
- **n = 13** (block 64): symplectic; **unique** — solved on the
  $+$ block: symmetric dimension 0, antisymmetric dimension 1, at two
  primes (the two blocks are exchanged by any single $\gamma_a$, which
  transports invariant forms between them).
- **n = 14** (blocks 64+64): **no invariant bilinear form on a single
  block, of either type** — solved on **both** blocks: symmetric 0,
  antisymmetric 0, at two primes. The two blocks are dual to each other
  through the full-module symmetric form (which swaps them). Same shape
  as $n = 10$.
- **n = 15** (block 128): symmetric (restriction of $J$); uniqueness via
  the same Schur argument as n = 13.
- **n = 16** (blocks 128+128): symmetric on each block (restriction of
  $J$). Uniqueness: direct elimination solve attempted three times and
  OOM-killed each time (see §4); uniqueness therefore rests on the Schur
  argument, labeled as such, not on a completed solve.

## 3. How the spinor checks were done (exactness and validation)

Explicit dense Sage matrices at module dimensions 128/256 with up to
7,140 homomorphism pairs over QQbar/ℚ(i) were judged prohibitively slow,
so the primary spinor engine represents every Clifford monomial **exactly**
as a Pauli string $(f, c, \mu)$ — flip bitmask, phase in $\mathbb Z/4$,
phase-pattern in $(\mathbb Z/4)^m$ — with the operator law
$e_x \mapsto i^{\,c + \sum_j \mu_j \mathrm{bit}_j(x)} e_{x \oplus f}$.
Multiplication, trace, transpose and adjoint are exact integer arithmetic
(no floating point, no approximation). Scripts: `pauli_monomial.py`,
`so11_so16_spinor.py`.

The engine is tied to explicit matrices by an independent validation web:

- **Control n = 10 (dense Sage, ℚ(i))**: gamma anticommutation on all
  ordered pairs: True; homomorphism on all 990 basis pairs: 0 failures;
  span rank 45; split 16+16 — the campaign's so(10) rung reproduced by
  dense matrices built from the monomial data. A smoke test also verified
  $mdense(\gamma_a\gamma_b) = mdense(\gamma_a)\,mdense(\gamma_b)$.
- **Dense Sage homomorphism, all pairs, exact**: n = 11 (1485 pairs,
  0 failures), n = 12 (2145 pairs, 0 failures).
  Script: `so11_so16_spinor_dense.py` (modes `hom`, `control10`).
- **Dense Sage exact samples**: n = 13, 14, 15, 16 — 50 randomly sampled
  pairs each (seed 20261016 + n), 0 failures (mode `sample`).
- **Numerical all-pairs (numpy, complex128)**: n = 14 (4095 pairs),
  n = 15 (5460), n = 16 (7140) — max residual **exactly 0.0** at every
  rung: the dense matrices have entries in $\{0, \pm1, \pm i\}/2$, exactly
  representable, so this pass is exact in effect; it is labeled CHECKED
  (numeric route) regardless. Script: `so11_so16_spinor_numeric.py`.
- **Span rank certificate (Sage)**: rank of $\{2\Sigma_{ab}\}$ over
  GF(p) at **two** split primes $p = 1000033$ and $p = 104729$ equals
  $\dim \mathfrak{so}(n)$ at every rung. Since the entries lie in
  $\mathbb Z[i]$, a maximal nonzero minor mod a prime ideal above $p$ is
  a nonzero minor over $\mathbb Z[i]$; hence the rank over ℚ(i) is
  maximal — the span is PROVED, not sampled (mode `ranks`; the same run
  re-verified the Γ split dimensions via kernel ranks mod p, and the
  vector Casimirs over QQ).
- **Form uniqueness solves (Sage, GF(p), two primes)**: on a Γ-eigenspace
  block, parametrize all symmetric (resp. antisymmetric) matrices $X$ and
  solve $\rho(e_{a,a+1})^T X + X \rho(e_{a,a+1}) = 0$ for the $n-1$
  adjacent generators (a generating set; invariance under it implies
  invariance under the whole algebra). Solution dimensions are as stated
  in §2. Scripts: `so11_so16_form_unique.py` (n = 11, 12, direct
  assembly), `so11_so16_form_unique2.py` (iterative-kernel assembly;
  n = 12 − block, n = 13 + block, n = 14 both blocks).

Status labels, rung by rung: path closures PROVED (exact) at all rungs;
spinor homomorphism PROVED by exact engine arithmetic at all rungs and
independently PROVED by dense Sage at n = 10, 11, 12 (all pairs) with
exact dense samples at 13–16 and numeric all-pairs (residual 0) at 14–16;
span ranks PROVED (two-prime certificate); splits and Casimir values
PROVED (exact); form existence/type PROVED (exact intertwiner identities
on all pairs); form uniqueness PROVED by solve at n = 11, 12, 13 and
nonexistence PROVED by solve at n = 14 — uniqueness at n = 15, 16 is
the standard Schur argument (labeled, not computed).

## 4. Disclosures — failures, fallbacks, environment incidents

1. **First path attempt lost.** The path closures were first launched as
   one laptop-side background job; the WSL distro restarted mid-run and
   wiped `/tmp` before any rung finished (no results captured). Rerun
   rung-by-rung with logs captured on this VM; all six rungs then
   completed (§1). No other work was affected.
2. **Two Sage-stdin bugs found and fixed.** (a) Comment lines inside
   function bodies split piped Sage blocks (`SyntaxError: incomplete
   input`); piped scripts were stripped of comment/blank lines.
   (b) Sage's preparser rewrites `^` as `**`, so `x ^ f` (xor) in the
   dense-matrix builder became a power and crashed with `IndexError`;
   replaced by an arithmetic `bxor`. The fixed builder passed the smoke
   test in §3 before the full runs.
3. **n = 16 form-uniqueness solve OOM-killed three times.** The block is
   128-dimensional (8,256 symmetric unknowns); the dense elimination
   process reached ~3.3 GB RSS and was killed by the OOM killer all three
   times (verified in `dmesg`: "Out of memory: Killed process …
   (sage-ipython)"), including once with the laptop otherwise idle. It
   was not re-attempted with a sparse/streaming variant; §2 labels the
   n = 16 uniqueness accordingly. All other computations completed.
4. **Fallback ladder.** The task's protocol asked for exact Sage
   homomorphism all-pairs at n = 11–13 and allowed numeric+sample at
   14–16. Delivered: exact dense Sage all-pairs at 11, 12 (13 covered by
   the exact engine all-pairs + dense exact 50-pair sample + two-prime
   rank certificate); at 14–16, exact engine all-pairs + dense exact
   samples + numeric all-pairs with residual exactly 0. No rung rests on
   numerics alone.
