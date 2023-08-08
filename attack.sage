load("sidh.sage")
load("uvtable.sage")
load("richelot_isogeny.sage")

class CDParams:
    """
    - start_sidh_pub: 攻撃対象の SIDH の公開情報
    - E: Bob の公開鍵 E_b
    - P, Q: φ(Pa), φ(Qa) (φ は Bob の秘密鍵)
    
    - iter_prime: 攻撃イテレーション i での 2^(a - α_(i-1)), 3^(b - β_(i - 1)) の情報
    - proof: パラメータのチェック
    """
    def __init__(self, start_sidh_pub: SIDHPublic, E, P, Q, iter_prime=None, betas=[0], ks=[], proof=True):
        if iter_prime == None:
            iter_prime = start_sidh_pub.prime

        E_start = start_sidh_pub.E0

        if proof:
            assert E.base_ring() == E_start.base_ring()
            assert P in E
            assert Q in E
            assert P.order() == 2^iter_prime.a
            assert Q.order() == 2^iter_prime.a
            assert P.weil_pairing(Q, 2^iter_prime.a) != 1

            assert E_start.a_invariants() in [(0, 0, 0, 1, 0), (0, 6, 0, 1, 0)]

            assert len(betas) == len(ks) + 1
        
        self.prime = prime
        self.start_sidh_pub = start_sidh_pub
        self.E_start = E_start

        self.E = E
        self.P = P
        self.Q = Q

        self.iter_prime = iter_prime
        self.betas = betas
        self.ks = ks

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

# i_isogeny(i_isogeny(P)) == -P
def i_isogeny(sidh_pub, P):
    E_start = sidh_pub.E0
    prime = sidh_pub.prime
    assert P in E_start

    Fpp = E_start.base_ring()
    x, y = P.xy()
    Q = E_start(-x, i*y)
    return Q

def is_split(C, E, Pc, P, Qc, Q, ai):
    h, D2_PcP, D2_QcQ = FromProdToJac(C, E, Pc, P, Qc, Q, ai)

    hp, D1, D2 = h, D2_PcP, D2_QcQ
    for i in range(1, ai - 1):
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

# 3^iter_prime.b isogeny φ:E1 -> E で，φ(P1) = P，φ(Q1) = Q となるものは存在するか？
def solve_D(params: CDParams, ui, vi, ker_kappa_gen, E1, P1, Q1):
    iter_prime = params.iter_prime
    ai = iter_prime.a
    bi = iter_prime.b
    alpha = params.start_sidh_pub.prime.a - ai
    beta = params.start_sidh_pub.prime.b - bi
    assert 2^ai - 3^bi == ui^2 + 4*vi^2
    assert ker_kappa_gen in params.E_start

    # calculate aixiliary isogeny
    assert params.E_start.a_invariants() == (0, 0, 0, 1, 0), "Not implemented (0, 6, 0, 1, 0) auxiliary isogeny"

    def gamma_start(P):
        assert P in params.E_start
        iP = i_isogeny(params.start_sidh_pub, P)
        assert i_isogeny(params.start_sidh_pub, iP) == -P
        return ui * P + (2 * vi) * iP

    ker_tilde_kappa_hat_gen = gamma_start(ker_kappa_gen)
    tilde_kappa_hat = params.E_start.isogeny(ker_tilde_kappa_hat_gen, algorithm="factored")
    C = tilde_kappa_hat.codomain()
    
    def gamma(P):
        assert P in params.E_start
        Q = 2^alpha * tilde_kappa_hat(gamma_start(P))
        return Q
    
    Pc = gamma(params.start_sidh_pub.Pa)
    Qc = gamma(params.start_sidh_pub.Qa)

    return is_split(C, params.E, Pc, params.P, Qc, params.Q, ai)

def attack(params: CDParams, iteration=1):
    assert len(params.betas) == iteration

    prime = params.prime
    iter_prime = params.iter_prime
    prev_ai = iter_prime.a
    prev_bi = iter_prime.b

    # last brute force iteration
    if prev_bi <= 3:
        print(f"[Last Iteration]")
        print(f"brute force 3^{prev_bi} kappa")
        print(f"betas: {params.betas}")
        print(f"ks: {params.ks}")
        assert len(params.betas) == len(params.ks) + 1
        recovered_skb = 0
        for i in range(len(params.ks)):
            recovered_skb += params.ks[i] * 3^params.betas[i]
        return recovered_skb

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
            assert 2^ai - 3^bi == ui^2 + 4*vi^2
            break
        alpha = prime.a - ai
        beta = prime.b - bi

    print(f"a_{iteration}: {ai}")
    print(f"b_{iteration}: {bi}")
    print(f"alpha_{iteration}: {alpha}")
    print(f"beta_{iteration}: {beta}")
    print()
    next_iter_prime = SIDHPrime(ai, bi, iter_prime.f, proof=False)

    # ki = 0, 1, 2 とかだと zero division error がおきる．なぜ
    ki = 0
    while True:
        kappa, ker_kappa_gen = choice_kappa(params, beta, ki)
        E1 = kappa.codomain()
        assert kappa.domain() == params.E_start
        assert kappa.degree() == 3^beta
        print(kappa)
        print(E1)
        print()

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
            break
        ki += 1

    return attack(next_params, iteration=iteration+1)
