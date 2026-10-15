# Wick rotation via the sqrt(FlatWave) complex structure J — a negative result

Date: 2026-10-15. Question (Kit): the pair of
`sqrt_fw_clifford.md` §2 supplies a canonical complex structure
$J = SG$ with $J^2 = -M_{FW}$ on the FlatWave $= +1$ function
sectors. Can that $J$ perform the Wick rotation of the
$\mathfrak{sl}(3,\mathbb{R})$ family — i.e. replace the imported
$i$ of E101 — when applied to the *full* basis (not a generating
set, per the Wick wall of `wick_sl8_wall.md`)?

Answer: **NO** (Sage-exact over QQbar; script
`proofs/wick_via_j.py`).

## Test design

Family: the E96 samples at the standard 10 degrees (QI + QII),
which close to $\mathfrak{sl}(3,\mathbb{R})$ (dim 8, re-verified
in the same run).

**Test arm.** Adjoin the $J$-partners of the samples. For a
$\pi$-periodic matrix function $F$, $JF = -g \cdot F$ pointwise
($g$ = the accumulated root, antiperiodic over $\pi$; $F(x+\pi) =
F(x)$). On $FW = +1$ sectors $g = \pm 1$, so $J$ collapses to a
sign there and contributes nothing; the only new generators are
$i \times$ (the QII samples). Compute the real Lie closure
inside $\mathfrak{sl}(3,\mathbb{C})$.

**Control arm.** Take the closed 8-dim basis, split it by matrix
transpose (3 antisymmetric + 5 symmetric, verified), and apply
the standard Wick rotation ($i$ on the symmetric five), as in
E101. Compute the real closure and the Killing signature.

## Results

| Arm | Real closure dim | Killing |
|---|---|---|
| Real samples alone | 8 | (5,3,0) — $\mathfrak{sl}(3,\mathbb{R})$ |
| Test: samples + $J$-partners | **16** | — (not 8-dim) |
| Control: standard Wick on the basis | 8 | **(0, 8, 0)** — $\mathfrak{su}(3)$ |

The test arm's 16 is the real dimension of
$\mathfrak{sl}(3,\mathbb{C})$ itself: the $J$-augmented family
generates the whole complexified algebra as a real algebra, not
any 8-dimensional real form. The control reproduces E101's
$\mathfrak{su}(3)$ in the same framework, so the apparatus does
detect (0,8) when a construction produces it.

## Verdict

$J$ does not perform the Wick rotation of the
$\mathfrak{sl}(3)$ family. The obstruction is categorical:
$J$ is a complex structure on the *function* space — it acts on
the argument $x$ (quarter/half-turn transport and multiplication
by $g$) — while the Wick rotation acts on the *value* space,
splitting the matrices by transpose, about which $J$ carries no
information. What $J$ achieves is the complexification
$\mathfrak{sl}(3,\mathbb{R}) \to
\mathfrak{sl}(3,\mathbb{C})_{\mathbb{R}}$; selecting a compact
real *form* inside the complexification requires a conjugation
(Cartan involution) in addition to a complex structure — the
same two-real-structures phenomenon recorded for the 1820
$\varepsilon$-refinement (Volume I audit workspace,
`io_word_sweep/EPS_ALIGNMENT_GAMMA_ROOT.md` §5).

WALL-CPLX bookkeeping: $J$'s canonical existence on the
$FW = +1$ function sectors stands (sqrt_fw_clifford §2, the
2026-10-15 Wick examination); the route "Wick-rotate the
$\mathfrak{sl}(3)$ family via $J$" is FALSIFIED. E101's Wick
rotation continues to rest on the imported $i$ and the
transpose split.

Scope: exact computation for the specific $J$ of §2 and the
E96 family. No statement is made about other complex
structures.
