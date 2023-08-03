load("sidh.sage")
load("uvtable.sage")

class CDParams:
    """
    - start_sidh_pub: 攻撃対象の SIDH の公開情報
    - E: Bob の公開鍵 E_b
    - P, Q: φ(Pa), φ(Qa) (φ は Bob の秘密鍵)
    
    - iter_prime: 攻撃イテレーション i での 2^(a - α_(i-1)), 3^(b - β_(i - 1)) の情報
    - proof: パラメータのチェック
    """
    def __init__(self, start_sidh_pub: SIDHPublic, E, P, Q, iter_prime=None, proof=True):
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
        
        self.start_sidh_pub = start_sidh_pub
        self.E_start = E_start

        self.E = E
        self.P = P
        self.Q = Q

        self.iter_prime = iter_prime
        self.prime = prime

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

def attack(params: CDParams, iteration=1):
    prime = params.prime
    iter_prime = params.iter_prime
    prev_ai = iter_prime.a
    prev_bi = iter_prime.b

    if prev_bi <= 3:
        print(f"[Last Iteration]")
        print(f"brute force 3^{prev_bi} kappa")
        return

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
    next_params = CDParams(params.start_sidh_pub, params.E, params.P, params.Q, iter_prime=next_iter_prime, proof=False)
    attack(next_params, iteration=iteration+1)
