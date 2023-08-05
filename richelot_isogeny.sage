# From A Note on Reimplementing the Castryck-Decru Attack and Lessons Learned for SageMath
def FromProdToJac(C, E, Pc, P, Qc, Q, ai):
    K = E.base_ring()
    assert C.base_ring() == K
    assert K.characteristic() != 2
    assert Pc in C
    assert P in E
    assert Qc in C
    assert Q in E
    assert Pc.order() == 2^ai
    assert P.order() == 2^ai
    assert Qc.order() == 2^ai
    assert Q.order() == 2^ai

    Pc2 = 2^(ai - 1) * Pc
    P2 = 2^(ai - 1) * P
    Qc2 = 2^(ai - 1) * Qc
    Q2 = 2^(ai - 1) * Q

    PR.<x> = PolynomialRing(K)

    alphas = [Pc2[0], Qc2[0], (Pc2 + Qc2)[0]]
    betas = [P2[0], Q2[0], (P2 + Q2)[0]]

    M = [
            [alphas[0] * betas[0], alphas[0], betas[0]],
            [alphas[1] * betas[1], alphas[1], betas[1]],
            [alphas[2] * betas[2], alphas[2], betas[2]]
        ]
    M = Matrix(K, M)
    D = M.determinant()
    [[R], [S], [T]] = M.inverse() * Matrix(K, [[1], [1], [1]])

    for i in range(3):
        assert (R * alphas[i] + T) * (R * betas[i] + S) == R + S * T

    delta_alpha = (alphas[0] - alphas[1]) * (alphas[1] - alphas[2]) * (alphas[2] - alphas[0])
    delta_beta = (betas[0] - betas[1]) * (betas[1] - betas[2]) * (betas[2] - betas[0])
    ss = [-delta_alpha / (R * D), -T / R]
    ts = [delta_beta / (R * D), -S / R]
    assert R + S * T == R^2 * ss[0] * ts[0]

    h_alphas = []
    h = ss[0]
    for i in range(3):
        h_alpha = (alphas[i] - ss[1]) / ss[0]
        h_alphas.append(h_alpha)
        h *= (x^2 - h_alpha)
    H = HyperellipticCurve(h)
    JH = H.jacobian()

    def Phi_hat(Pc, Pe):
        print("yey")
        assert Pc in C
        assert Pc.order() == 2^ai
        assert Pe in E
        assert Pe.order() == 2^ai

        # phi1_hat: C -> JH
        x1, y1 = Pc.xy()
        U = ss[0] * x^2 + ss[1] - x1
        V = PR(y1 / ss[0])
        JPc = JH([U, V])

        # phi2_hat: E -> JH
        x2, y2 = Pe.xy()
        U = ts[0] - x^2 * (x2 - ts[1])
        V = x^3 * (y2 / ts[0])
        JPe = JH([U, V])
        return JPc + JPe

    D2_PcP = Phi_hat(Pc, P)
    D2_QcQ = Phi_hat(Qc, Q)

