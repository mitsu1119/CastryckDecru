def FromProdToJac(C, E, Pc, P, Qc, Q, ai):
    K = E.base_ring()
    assert C.base_ring() == K
    assert K.characteristic() != 2

    PR.<x> = PolynomialRing(K)

    # factoring C/K: y^2 = (x - alphas[0])(x - alphas[1])(x - alphas[2])
    coeffs = C.a_invariants()
    assert (coeffs[0], coeffs[2]) == (0, 0)
    FC = x^3 + coeffs[1] * x^2 + coeffs[3] * x + coeffs[4]
    alphas = []
    for root, e in FC.roots():
        alphas += [root] * e
    assert len(alphas) == 3
    assert FC == (x - alphas[0]) * (x - alphas[1]) * (x - alphas[2])

    # factoring E/K: y^2 = (x - betas[0])(x - betas[1])(x - betas[2])
    coeffs = E.a_invariants()
    assert (coeffs[0], coeffs[2]) == (0, 0)
    FE = x^3 + coeffs[1] * x^2 + coeffs[3] * x + coeffs[4]
    betas = []
    for root, e in FE.roots():
        betas += [root] * e
    assert len(betas) == 3
    assert FE == (x - betas[0]) * (x - betas[1]) * (x - betas[2])

