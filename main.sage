load("attack.sage")

a = 33
b = 19
f = 1
prime = SIDHPrime(a, b, f)

Fpp.<i> = GF(prime.p^2, modulus=x^2+1)
E0 = EllipticCurve(Fpp, [0, 0, 0, 1, 0])

Pa, Qa, Pb, Qb = generate_sidh_torsions(E0, prime)
sidh_pub = SIDHPublic(prime, E0, Pa, Qa, Pb, Qb)
skb = 3271984327985
phi_b = E0.isogeny(Pb + skb * Qb, algorithm="factored")

print("[target secrets]")
print(f"secret: skb = {skb}")
print(phi_b)
print()

E = phi_b.codomain()
P = phi_b(Pa)
Q = phi_b(Qa)

bobs_public = BobsPublic(sidh_pub, E, P, Q)
cdparams = CDParams(sidh_pub, E, P, Q)
recovered_skb = attack(cdparams, bobs_public)

print()
assert recovered_skb == skb
print(f"recovered!! skb = {recovered_skb}")
