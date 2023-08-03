class SIDHPrime:
    def __init__(self, a, b, f, proof=True):
        self.a = a
        self.b = b
        self.f = f
        self.p = 2^a * 3^b * f - 1

        if proof:
            assert is_prime(self.p)
            assert gcd(6, f) == 1

    def __str__(self):
        return f"SIDHPrime {self.p} = 2^{self.a} * 3^{self.b} * {self.f} - 1"
