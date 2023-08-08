class SIDHPrime:
    def __init__(self, a, b, f, proof=True):
        p = 2^a * 3^b * f - 1

        if proof:
            assert is_prime(p)
            assert gcd(6, f) == 1

        self.a = a
        self.b = b
        self.f = f
        self.p = p

    def __str__(self):
        return f"SIDHPrime {self.p} = 2^{self.a} * 3^{self.b} * {self.f} - 1"

class SIDHPublic:
    def __init__(self, prime: SIDHPrime, E0, Pa, Qa, Pb, Qb, proof=True):
        if proof:
            Fpp.<i> = GF(prime.p^2, modulus=x^2+1)
            assert E0.base_ring() == Fpp
            assert E0.order() == (prime.p + 1)^2

            assert Pa in E0
            assert Qa in E0
            assert Pa.order() == 2^prime.a
            assert Qa.order() == 2^prime.a
            assert Pa.weil_pairing(Qa, 2^prime.a) != 1

            assert Pb in E0
            assert Qb in E0
            assert Pb.order() == 3^prime.b
            assert Qb.order() == 3^prime.b
            assert Pb.weil_pairing(Qb, 3^prime.b) != 1

        self.prime = prime
        self.E0 = E0
        self.Pa = Pa
        self.Qa = Qa
        self.Pb = Pb
        self.Qb = Qb

class BobsPublic:
    def __init__(self, sidh: SIDHPublic, Eb, P, Q):
        assert Eb.base_ring() == sidh.E0.base_ring()
        assert P in Eb
        assert Q in Eb
        assert P.order() == 2^sidh.prime.a
        assert Q.order() == 2^sidh.prime.a
        assert P.weil_pairing(Q, 2^sidh.prime.a) != 1

        self.sidh = sidh
        self.Eb = Eb
        self.P = P
        self.Q = Q

def generate_sidh_torsions(E0, prime: SIDHPrime):
    Fpp.<i> = GF(prime.p^2, modulus=x^2+1)
    assert E0.base_ring() == Fpp
    assert E0.order() == (prime.p + 1)^2

    k = (prime.p + 1) // 2^prime.a
    Pa, Qa = [k * G for G in E0.gens()]
    k = (prime.p + 1) // 3^prime.b
    Pb, Qb = [k * G for G in E0.gens()]

    assert Pa.order() == 2^prime.a
    assert Qa.order() == 2^prime.a
    assert Pb.order() == 3^prime.b
    assert Qb.order() == 3^prime.b
    assert Pa.weil_pairing(Qa, 2^prime.a) != 1
    assert Pb.weil_pairing(Qb, 3^prime.b) != 1
    return Pa, Qa, Pb, Qb
