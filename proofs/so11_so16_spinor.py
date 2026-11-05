import sys, time
sys.path.insert(0, '.')
from pauli_monomial import (mmul, mneg, mscale_i, ident, is_scalar, trace_is_zero,
                      mtranspose, madjoint, gamma, mprod_list, so_bracket)

def analyse(n, control=False):
    t0 = time.time()
    even = (n % 2 == 0)
    m = n // 2 if even else (n + 1) // 2
    d = 2 ** m
    I = ident(m)
    gam = dict((a, gamma(a, m)) for a in range(1, n + 1))
    # anticommutation + hermiticity of gammas
    anti_ok = True
    for a in range(1, n + 1):
        if madjoint(gam[a]) != gam[a]:
            anti_ok = False
        if mmul(gam[a], gam[a]) != I:
            anti_ok = False
        for b in range(a + 1, n + 1):
            if mmul(gam[a], gam[b]) != mneg(mmul(gam[b], gam[a])):
                anti_ok = False
    pairs = [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]
    M = dict((p, mmul(gam[p[0]], gam[p[1]])) for p in pairs)
    # Sigma^dagger = -Sigma and Sigma^2 = -I/4 (Casimir input)
    cas_ok = all(madjoint(M[p]) == mneg(M[p]) and mmul(M[p], M[p]) == (0, 2, (0,) * m)
                 for p in pairs)
    # homomorphism, engine-exact, all unordered distinct basis pairs
    hom_fail = []
    npairs = len(pairs)
    for ii in range(npairs):
        for jj in range(ii + 1, npairs):
            (a, b), (c, dd) = pairs[ii], pairs[jj]
            P = mmul(M[(a, b)], M[(c, dd)])
            Q = mmul(M[(c, dd)], M[(a, b)])
            br = so_bracket(a, b, c, dd)
            if not br:
                ok = (P == Q)
            else:
                assert len(br) == 1
                (xy, sgn), = br.items()
                ok = (P == mneg(Q)) and (P == (M[xy] if sgn == 1 else mneg(M[xy])))
            if not ok:
                hom_fail.append(((a, b), (c, dd)))
    # chirality / volume element
    w = mprod_list([gam[a] for a in range(1, n + 1)])
    w2 = mmul(w, w)
    assert is_scalar(w2) and all(v % 4 == 0 for v in w2[2])
    G = w if w2[1] % 4 == 0 else mscale_i(w)
    G_ok = (mmul(G, G) == I) and trace_is_zero(G) and \
        all(mmul(M[p], G) == mmul(G, M[p]) for p in pairs)
    # rank: Gram matrices (full module and G-blocks) diagonal?
    def gram_offdiag_zero(withG):
        for ii in range(npairs):
            for jj in range(ii + 1, npairs):
                prod = mmul(M[pairs[ii]], M[pairs[jj]])
                if withG:
                    prod = mmul(G, prod)
                if is_scalar(prod):
                    return False
        return True
    rank_full = npairs if gram_offdiag_zero(False) else -1
    rank_block = npairs if (gram_offdiag_zero(False) and gram_offdiag_zero(True)) else -1
    # invariant bilinear form: monomial charge-conjugation search
    evens = [a for a in range(1, n + 1) if a % 2 == 0]
    odds = [a for a in range(1, n + 1) if a % 2 == 1]
    cand = [evens, odds]
    if not even:
        cand = cand + [evens + [n + 1], odds + [n + 1]]
    gamfull = dict(gam)
    if not even:
        gamfull[n + 1] = gamma(n + 1, m)
    eps = dict((cc, 1 if cc % 2 == 1 else -1) for cc in range(1, n + 1))
    form_info = None
    for S in cand:
        J = mprod_list([gamfull[a] for a in S])
        def tau_test(JJ):
            ss = []
            for cc in range(1, n + 1):
                lhs = mmul(JJ, gam[cc])
                rhs = mmul(gam[cc], JJ)
                if lhs == rhs:
                    ss.append(1)
                elif lhs == mneg(rhs):
                    ss.append(-1)
                else:
                    return None
            for sigma in (1, -1):
                if all(ss[k] == sigma * eps[k + 1] for k in range(n)):
                    return sigma
            return None
        sigma = tau_test(J)
        if sigma is None:
            continue
        JG, GJ = mmul(J, G), mmul(G, J)
        if JG == GJ:
            block = 'preserves-blocks'
        elif JG == mneg(GJ):
            block = 'swaps-blocks'
            if not even:
                # multiply by the volume element: preserves blocks, same tau type
                J = mmul(J, G)
                if tau_test(J) is None:
                    continue
                JG, GJ = mmul(J, G), mmul(G, J)
                if JG != GJ:
                    continue
                block = 'preserves-blocks (after x Gamma_vol)'
        else:
            continue
        Jt = mtranspose(J)
        if Jt == J:
            ftype = 'symmetric'
        elif Jt == mneg(J):
            ftype = 'antisymmetric (symplectic)'
        else:
            ftype = 'NEITHER'
        full_ok = all(mmul(mtranspose(M[p]), J) == mneg(mmul(J, M[p])) for p in pairs)
        form_info = (S, sigma, ftype, block, full_ok)
        break
    print('SPIN n=%d m=%d moddim=%d anti=%s casimirM2=%s hom_fail=%d rank=%d blockrank=%d '
          'G_ok=%s split=%d+%d form=%s elapsed=%.1fs'
          % (n, m, d, anti_ok, cas_ok, len(hom_fail), rank_full, rank_block,
             G_ok, d // 2, d // 2, form_info, time.time() - t0), flush=True)
    if hom_fail:
        print('  HOM FAILURES:', hom_fail[:5], flush=True)

def main():
    for n in [10, 11, 12, 13, 14, 15, 16]:
        analyse(n)
    print('SPIN ENGINE ALL DONE', flush=True)

main()
