# The Spinor Bundle over OG+(5,10), Built

Date: 2026-10-23. Status: PROVED (explicit exact matrices over Q
and Q(i), 20 counted checks, exit 0). Script:
proofs/spinor_bundle.py (sympy, this machine; ~3 s). Third of the
16's layer builds, on pure_spinor_orbit.md (the base) and
su5_stabilizer.md (the isotropy group and branching, imported
here as proved input).

## What "the bundle" is, concretely

OG+(5,10) is the homogeneous space Spin(10)/U(5) (compact form:
SO(10)/U(5) up to the cover). The spinor bundle is the associated
bundle Spin(10) x_{U(5)} 16. A homogeneous bundle is determined
by its isotropy representation, so the buildable content is:
the tangent model of the base (the reductive complement and its
isotropy action), the spinor fiber's isotropy action (proved in
the stabilizer build), the tautological line bundle, and the
descent question -- whether the bundle is a Spin bundle or
secretly an SO(10) bundle. All four are settled here by exact
computation.

## Verified exactly (script sections)

1. Reductive complement. m = {real so(10) elements with adapted
   A-block = 0} is SOLVED (real nullspace): dim_R m = 20.
   h (+) m = so(10) (rank 45). Brackets on full bases:
   [h, m] subset m (25 x 20 pairs -- m is an isotropy module, the
   tangent space of the base at the basepoint) and
   [m, m] subset h (all 190 pairs) -- the Hermitian-symmetric
   relation; SO(10)/U(5) is the type-DIII space. dim_R 20 =
   2 * dim_C 10, matching the orbit build's complex dimension.
   Controls: the explicit B = C element is in m (true);
   E_{0,2} is not (false -- its A-block is nonzero).
2. Complexification. The B- and C-projections of m_C each span
   all of Lambda^2 (complex rank 10 each; the reality link
   C = conj(B) for real elements makes each projection a complex
   subspace, so the doubled-real-rank count is exact). Hence
   m_C = n (+) n-bar with n the B-block radical of the orbit
   build's parabolic: the (1,0)/(0,1) split of the tangent.
3. Tangent isotropy action. For every stabilizer basis element
   (adapted block diag(A, -A^T)) and every radical basis element
   (B running over a complex basis of skew matrices):
   ad sends B to A B + B A^T, all other blocks of the commutator
   zero (25 x 20 pairs). In B-matrix coordinates this is exactly
   the derivation action of gl(W) on Lambda^2: the holomorphic
   tangent at the basepoint is the Lambda^2 isotropy module.
4. Descent (Spin, not SO). The centre generator Z = -sum_k
   E_{2k,2k+1} is verified central in h, with A(Z) = -i I,
   D(Z) = +i I. rho(Z) is diagonal on the 32 with every eigenvalue
   in i*{+-1/2, +-3/2, +-5/2}. Hence exp(2 pi rho(Z)) = -I_32
   (each scalar exponential is e^{+-i pi * odd} = -1) while
   exp(2 pi Z) on vectors is +I_10 (e^{+-2 pi i} = 1 on both
   blocks). The circle closes at 2 pi in SO(10) and at 4 pi in
   Spin(10); the -I_32 is the covering kernel acting. The
   isotropy representation does not descend to SO(10)'s U(5):
   the spinor bundle is a genuine Spin bundle.
5. Tautological line bundle L (fiber over a point of OG+ = the
   pure spinor line itself). Its isotropy character is the
   mu = (1/2) tr A of the orbit/stabilizer builds (re-verified
   from the spinor column on all 25 basis elements); its Z-weight
   is -5/2, i.e. the singlet summand of the spinor bundle IS L:
   as isotropy modules, S = L_mu (x) (Lambda^0 (+) Lambda^2 (+)
   Lambda^4 of the creation-space bundle), with the branching of
   su5_stabilizer.md.

## Scope

This is the homogeneous (universal) model: OG+ with its
tautological Spin(10)-structure. A Spin(10)-structure on an
arbitrary 10-manifold, and its spinor bundle, is the standard
associated-bundle construction applied to a reduction of the
frame bundle; no such manifold is given or constructed here, and
none is needed for the layer as listed. Global existence of the
homogeneous bundle from (G, H, isotropy rep) is the standard
associated-bundle theorem, cited as input; everything about the
fibers, tangent, characters, and descent is computed. No FW/saw
operator content. Local file; not pushed.
