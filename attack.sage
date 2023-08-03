load("sidh.sage")

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
        prime = iter_prime

        if proof:
            assert E.base_ring() == E_start.base_ring()
            assert P in E
            assert Q in E
            assert P.order() == 2^prime.a
            assert Q.order() == 2^prime.a
            assert P.weil_pairing(Q, 2^prime.a) != 1

            assert E_start.a_invariants() in [(0, 0, 0, 1, 0), (0, 6, 0, 1, 0)]
        
        self.start_sidh_pub = start_sidh_pub
        self.E_start = E_start

        self.E = E
        self.P = P
        self.Q = Q

        self.prime = prime
