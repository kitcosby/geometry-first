# Tables Manifest

**Status:** The human-readable CSV tables (critical subsets, diagonals) are in this directory. The large binary `.npy`/`.npz` files (~532MB total) are not in the repo (GitHub file-size limits); they are available in the canonical `~/workspace/tables/` store. Every binary has a corresponding CSV here that losslessly rebuilds it or summarizes it.

**E8 tables: EXCLUDED.** The E8-specific tables (30380 generators/weights, E8 adjoint, simple-root CSVs) are not in this repository. See README.md.

## Binary tables (not in repo; available on request)

| File | Content | Status |
|------|---------|--------|
| `gammas.npy` | 16×256×256 SO(16) Clifford generators | NC |
| `Pi_plus.npy`, `Pi_minus.npy` | 256×256 chiral projectors | NC |
| `B.npy` | 256×256 charge conjugation | NC |
| `U_Splus.npy` | 256×128 S⁺ basis | NC |
| `M_so16.npy` | 120×128×128 so(16) generators | NC |
| `singlet.npy` | 16384-vector true singlet | NC |
| `fourform_Bmapped.npy` | 1820×8256 B-mapped four-forms | NC |
| `C.npz` | 8256×8256 Casimir (sparse) | NC |
| `P1.npz`, `P1820.npz`, `P6435.npz` | Casimir projectors | NC |
| `P54_135.npy` | 135×135 projector onto (54,1) | SC |
| `Pi0_grade_1820.npy` etc. | 1820×1820 grade projectors | NC |
| `Me_reference_1820.npy` | 1820×1820 reference mass op | NC |
| `M_eps_1820.npy` | 1820×1820 epsilon operator | SC |
| `Me_ioword_1820.npy` | 1820×1820 IO-word mass op | SC |
| `C1/C2/C3_*_1820.npy` | 1820×1820 IO-code operators | SC |
| `Pk_fourier_1820.npy` | 8×1820×1820 Fourier projectors | SC |
| `J120_so16_on_1820.npz` | 120 so(16) generators on 1820 | SC |
| `W2to4_120x1820x120_int8.npy` | Wedge bridges | SC |
| `so16_viol_norms.npy` etc. | SO(16) violation/block norms | NC/CHECKED |
| `Gamma10_intertwiners.npz` | SO(10) intertwiners | — |

**NC** = numerically checked. **SC** = symbolically checked. **CHECKED** = numerical verification.

## CSV tables (in this directory)

All `_critical.csv` and `_diag.csv` files: lossless critical subsets or diagonals. Each rebuilds or summarizes its binary exactly. See TABLES.md in the source catalog for provenance.
