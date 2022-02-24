#BLS6 curve taken from https://eprint.iacr.org/2019/431.pdf Table 15

x = -0xdffffffffffffffffbffffffffffffffffffffdfffffffc00001 
r = x^2 - x +1   
p = 1/3*(x-1)^2*(x^2 -x+1) +x 
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
D = fundamental_discriminant(-3)
H = FpX(hilbert_class_polynomial(D))
E = EllipticCurve_from_j(H.roots()[0][0])  
assert E.order() % r  == 0 

G1 = (E.order()//r) * E.random_point()
while G1.is_zero() :
    G1 = (E.order()//r) * E.random_point()

beta = Fp.random_element()
while not((X**6-beta).is_irreducible()):
    beta = Fp.random_element()
Fp6.<u> = GF(p**6, modulus = X**6-beta)
E6 = E.change_ring(Fp6)

i=1
E_t = EllipticCurve([E.a4(), E.a6() * beta**i])
while E_t.order()%r!= 0:
    i+=1
    E_t = EllipticCurve([E.a4(), E.a6() * beta**i])

G2_t = (E_t.order()//r)*E_t.random_point()
while G2_t.is_zero():
    G2_t = (E_t.order()//r)*E_t.random_point()

sqrt3_beta = u**2
sqrt_beta = u**3
G1 = E6(G1)
G2 = E6(G2_t[0] /sqrt3_beta, G2_t[1]/sqrt_beta)
eG1G2 = G1.tate_pairing(G2, r, 6)
