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
    FC = (x - alphas[0]) * (x - alphas[1]) * (x - alphas[2])
    FE = (x - betas[0]) * (x - betas[1]) * (x - betas[2])

    Delta_alpha = FC.discriminant()
    Delta_beta = FE.discriminant()

    a1 = (alphas[2] - alphas[1])^2 / (betas[2] - betas[1]) + (alphas[1] - alphas[0])^2 / (betas[1] - betas[0]) + (alphas[0] - alphas[2])^2 / (betas[0] - betas[1])
    b1 = (betas[2] - betas[1])^2 / (alphas[2] - alphas[1]) + (betas[1] - betas[0])^2 / (alphas[1] - alphas[0]) + (betas[0] - betas[2])^2 / (alphas[0] - alphas[1])
    a2 = alphas[0] * (betas[2] - betas[1]) + alphas[1] * (betas[0] - betas[2]) + alphas[2] * (betas[1] - betas[0])
    b2 = betas[0] * (alphas[2] - alphas[1]) + betas[1] * (alphas[0] - alphas[2]) + betas[2] * (alphas[1] - alphas[0])

    A = Delta_beta * a1 / a2
    B = Delta_alpha * b1 / b2

    h = -(A * (alphas[1] - alphas[0]) * (alphas[0] - alphas[2]) * x^2 + B * (betas[1] - betas[0]) * (betas[0] - betas[2])) * (A * (alphas[2] - alphas[1]) * (alphas[1] - alphas[0]) * x^2 + B * (betas[2] - betas[1]) * (betas[1] - betas[0])) * (A * (alphas[0] - alphas[2]) * (alphas[2] - alphas[1]) * x^2 + B * (betas[0] - betas[2]) * (betas[2] - betas[1]))

    H = HyperellipticCurve(h)
    JH = H.jacobian()

    s1 = -(B * a2) / (A * a1)
    s2 = (alphas[0] * (alphas[2] - alphas[1])^2 / (betas[2] - betas[1]) + alphas[1] * (alphas[0] - alphas[2])^2 / (betas[0] - betas[2]) + alphas[2] * (alphas[1] - alphas[0])^2 / (betas[1] - betas[0])) / a1
    t1 = -(A * b2) / (B * b1)
    t2 = (betas[0] * (betas[2] - betas[1])^2 / (alphas[2] - alphas[1]) + betas[1] * (betas[0] - betas[2])^2 / (alphas[0] - alphas[2]) + betas[2] * (betas[1] - betas[0])^2 / (alphas[1] - alphas[0])) / b1

    # Phi: JH -> CxE
    # Phi = phi_1.dual() x phi_2.dual()
    # find D in JH that is Phi(D) = (Pc, Pe) and calculate 2D
    # That mean this isogeny is (2, 2)-isogeny: CxE -> JH
    def Phi_hat(Pc, Pe):
        assert Pc in C
        assert Pe in E
        assert Pc.order() == 2^ai
        assert Pe.order() == 2^ai

        MumfordCoeffs.<U0, U1, V0, V1> = AffineSpace(Fpp, 4)

        # for phi1(PH + QH)
        U0_tilde = 1 / U0
        U1_tilde = U1 / U0
        V0_tilde = (U1 * V0 - U0 * V1) / U0^2
        V1_tilde = (U1^2 * V0 - U0 * V0 - U0 * U1 * V1) / U0^2
        lam1 = -(Delta_beta * V1_tilde) / (A^3 * s1 * U1_tilde)
        phi1x = lam1^2 + alphas[0] + alphas[1] + alphas[2] - s1 * (U1_tilde^2 - 2 * U0_tilde) - 2 * s2
        phi1y = -lam1 * (phi1x - s2 + (U0_tilde * V1_tilde - U1_tilde * V0_tilde) * (s1 / V1_tilde))

        # for phi2(PH + QH)
        lam2 = -(Delta_alpha * V1) / (B^3 * t1 * U1)
        phi2x = lam2^2 + betas[0] + betas[1] + betas[2] - t1 * (U1^2 - 2 * U0) - 2 * t2
        phi2y = -lam2 * (phi2x - t2 + (U0 * V1 - U1 * V0) * (t1 / V1))

        eq1 = phi1x - Pc[0]
        eq2 = phi1y - Pc[1]
        eq3 = phi2x - P[0]
        eq4 = phi2y - P[1]
        eq5 = 2*V0^2 - 2*V0*V1*U1 + V1^2*(U1^2 - 2*U0) - 2*h[0] - (-U1)*h[1] - (U1^2 - 2*U0)*h[2] - (-U1^3 + 3*U0*U1)*h[3] - (U1^4 - 4*U1^2*U0 + 2*U0^2)*h[4] - (-U1^5 + 5*U1^3*U0 - 5*U1*U0^2)*h[5] - (U1^6 - 6*U1^4*U0 + 9*U1^2*U0^2 - 2*U0^3)*h[6]

        eq1 = eq1.numerator()
        eq2 = eq2.numerator()
        eq3 = eq3.numerator()
        eq4 = eq4.numerator()
        eq5 = eq5.numerator()

        Id = ideal(eq1, eq2, eq3, eq4, eq5)

        # too late
        print(Id.groebner_basis())
        print(Id.variety())

        return None

    JPcP = Phi_hat(Pc, P)
