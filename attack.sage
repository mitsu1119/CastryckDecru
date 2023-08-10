from concurrent.futures.thread import ThreadPoolExecutor

load("CDParams.sage")
load("uvtable.sage")
load("richelot_isogeny.sage")

"""
prime := 2^a * 3^b * f - 1
c := 2^(a - alpha) - 3^(b - beta)
c が u^2 + 4*v^2 の形でかけるような最小の奇数 beta と，それに付随する alpha を探索
ただし b > 1
"""
def search_alpha_beta(prime: SIDHPrime):
    a = prime.a
    b = prime.b
    if b % 2 == 1:
        for beta in range((b % 2) + 1, b, 2):
            s = b - beta
            index = (s - 1) // 2
            assert s == uvtable[index][0]
            s, t, u, v = uvtable[index]
            if t <= a:
                alpha = a - t
                assert 2^(a - alpha) - 3^(b - beta) == u^2 + 4*v^2
                return alpha, beta, u, v
    return None

"""
復元したい同種写像 φ に含まれている可能性がある κ の候補を計算
choice によって κ の候補が変わる
"""
def choice_kappa(params: CDParams, beta, choice):
    prime = params.prime
    b = prime.b
    assert 1 <= beta <= b
    assert 0 <= choice <= 3^beta - 1

    ki = choice
    bi = b - beta

    assert len(params.ks) + 1 == len(params.betas)
    Q_coef = 0
    for i in range(len(params.ks)):
        Q_coef += params.ks[i] * 3^params.betas[i]
    Q_coef += ki * 3^params.betas[-1]

    R = 3^bi * (params.start_sidh_pub.Pb + Q_coef * params.start_sidh_pub.Qb)
    kappa = params.E_start.isogeny(R, algorithm="factored")

    return kappa, R

def is_split(C, E, Pc, P, Qc, Q, ai):
    h, D2_PcP, D2_QcQ = FromProdToJac(C, E, Pc, P, Qc, Q, ai)

    hp, D1, D2 = h, D2_PcP, D2_QcQ
    for i in range(1, ai - 1):
        print(f"{i}th Jac -> Jac")
        hp, D1, D2 = FromJacToJac(hp, D1, D2, ai - i)
        JHp = HyperellipticCurve(hp).jacobian()
        assert 2^(ai - i) * D1 == 0
        assert 2^(ai - i) * D2 == 0

    # last step
    res = FromJacToJac_last_test(hp, D1, D2, 1)
    
    # be split because delta == 0
    if res:
        return True

    # gluing
    return False

# isogeny 2i
# two_i_isogeny(two_i_isogeny(P)) == -4P
def two_i_isogeny(E, P):
    assert P in E
    if E.a_invariants() == (0, 0, 0, 1, 0):
        # i_isogeny(i_isogeny(P)) == -P
        def i_isogeny(PP):
            Fp2 = E.base_ring()
            x, y = PP.xy()
            Q = E(-x, i * y)
            return Q
        res = i_isogeny(P)
        assert i_isogeny(res) == -P
        return 2 * res

    if E.a_invariants() == (0, 6, 0, 1, 0):
        two_i_isogeny = E.isogeny(E.lift_x(ZZ(1)), codomain=E)
        res = two_i_isogeny(P)
        assert two_i_isogeny(res) == -4 * P
        return res

    assert 1 == 2, "Not implemented two_i_isogeny"

# 3^iter_prime.b isogeny φ:E1 -> E で，φ(P1) = P，φ(Q1) = Q となるものは存在するか？
def solve_D(params: CDParams, ui, vi, ker_kappa_gen, E1, P1, Q1):
    iter_prime = params.iter_prime
    ai = iter_prime.a
    bi = iter_prime.b
    alpha = params.start_sidh_pub.prime.a - ai
    beta = params.start_sidh_pub.prime.b - bi

    assert 2^ai - 3^bi == ui^2 + 4*vi^2
    assert ker_kappa_gen in params.E_start

    def gamma_start(P):
        assert P in params.E_start
        return ui * P + vi * two_i_isogeny(params.E_start, P)

    ker_tilde_kappa_hat_gen = gamma_start(ker_kappa_gen)
    tilde_kappa_hat = params.E_start.isogeny(ker_tilde_kappa_hat_gen, algorithm="factored")
    C = tilde_kappa_hat.codomain()
    
    # c-isogeny gamma: E_start -> C
    def gamma(P):
        assert P in params.E_start
        Q = 2^alpha * tilde_kappa_hat(gamma_start(P))
        return Q

    Pc = gamma(params.start_sidh_pub.Pa)
    Qc = gamma(params.start_sidh_pub.Qa)

    return is_split(C, params.E, Pc, params.P, Qc, params.Q, ai)

def push_correct_choiced_kappa(params: CDParams, prev_ai, next_iter_prime, ui, vi, alpha, beta, ki, kappa_buf):
    ai = next_iter_prime.a
    kappa, ker_kappa_gen = choice_kappa(params, beta, ki)
    E1 = kappa.codomain()
    assert kappa.domain() == params.E_start
    assert kappa.degree() == 3^beta

    P1 = kappa(2^alpha * params.start_sidh_pub.Pa)
    Q1 = kappa(2^alpha * params.start_sidh_pub.Qa)
    assert P1 in E1
    assert Q1 in E1

    assert P1.order() == 2^ai
    assert Q1.order() == 2^ai

    kappa_hat = kappa.dual()
    assert kappa_hat.domain() == E1

    alpha_diff = prev_ai - ai
    P_dest = 2^alpha_diff * params.P
    Q_dest = 2^alpha_diff * params.Q

    assert P_dest in params.E
    assert Q_dest in params.E
    assert P_dest.order() == 2^ai
    assert Q_dest.order() == 2^ai
    assert P_dest.weil_pairing(Q_dest, 2^ai) != 1

    next_params = CDParams(params.start_sidh_pub, params.E, P_dest, Q_dest, iter_prime=next_iter_prime, betas=params.betas+[beta], ks=params.ks+[ki])

    print(f"test ki = {ki}")
    if solve_D(next_params, ui, vi, ker_kappa_gen, E1, P1, Q1):
        print(f"skb = {ki} * 3^{beta}")
        print("split!!")
        print(f"ki: {next_params.ks}, betas: {next_params.betas}")
        print()
        kappa_buf.append(next_params)

def attack(params: CDParams, bobs_public: BobsPublic, iteration=1):
    assert len(params.betas) == iteration

    prime = params.prime
    iter_prime = params.iter_prime
    prev_ai = iter_prime.a
    prev_bi = iter_prime.b

    # last step
    if prev_bi <= 3:
        print(f"[Last Iteration]")
        print(f"brute force 3^{prev_bi} kappa")
        print(f"betas: {params.betas}")
        print(f"ks: {params.ks}")
        assert len(params.betas) == len(params.ks) + 1

        # last brute force iteration
        recovered_skb = 0
        for i in range(len(params.ks)):
            recovered_skb += params.ks[i] * 3^params.betas[i]

        last_beta = params.betas[-1]
        for ki in range(3^last_beta):
            skb_est = recovered_skb + ki * 3^last_beta
            R = params.start_sidh_pub.Pb + skb_est * params.start_sidh_pub.Qb
            phiB_est = params.start_sidh_pub.E0.isogeny(R, algorithm="factored")
            phiB_Pa = phiB_est(params.start_sidh_pub.Pa)
            phiB_Qa = phiB_est(params.start_sidh_pub.Qa)
            if phiB_Pa == bobs_public.P and phiB_Qa == bobs_public.Q:
                return skb
        assert 1 == 2, "brute force missing"

    # print the information for iteration
    if iteration == 1:
        print("[First Iteration]")
    else:
        print(f"[Iteration {iteration}]")
    print(f"target: 2^{prev_ai}, 3^{prev_bi}")

    # deciding alpha_i, beta_i, ai, bi 
    if iteration == 1:
        alpha, beta, ui, vi = search_alpha_beta(iter_prime)
        ai = prime.a - alpha
        bi = prime.b - beta
        assert 2^ai - 3^bi == ui^2 + 4*vi^2
    else:
        bi = prev_bi
        while True:
            bi -= 1
            table = uvtable_t(bi)
            if table == None:
                continue
            bi, ai, ui, vi = table
            if ai > prev_ai:
                continue
            assert 2^ai - 3^bi == ui^2 + 4*vi^2
            assert prev_ai >= ai
            break
        alpha = prime.a - ai
        beta = prime.b - bi

    print(f"a_{iteration}: {ai}")
    print(f"b_{iteration}: {bi}")
    print(f"alpha_{iteration}: {alpha}")
    print(f"beta_{iteration}: {beta}")
    print()
    next_iter_prime = SIDHPrime(ai, bi, iter_prime.f, proof=False)

    ki = 0
    kappa_next_params = []
    while True:
        push_correct_choiced_kappa(params, prev_ai, next_iter_prime, ui, vi, alpha, beta, ki, kappa_next_params)
        if len(kappa_next_params) >= 1:
            break
        ki += 1

    assert len(kappa_next_params) != 0
    next_params = kappa_next_params[0]
    return attack(next_params, bobs_public, iteration=iteration+1)
