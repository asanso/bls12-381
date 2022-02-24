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

Fp6.<u> = GF(p**6)
E6 = E.change_ring(Fp6)

G1 = (E.order()//r) * E.random_point()
while G1.is_zero() :
    G1 = (E.order()//r) * E.random_point()
G1 = E6(G1)

beta = Fp.random_element()
while len((X**6 - beta).roots())==0 :
    beta = Fp.random_element()

i=1
E_t = EllipticCurve([E.a4(), E.a6() * beta**i])
while E_t.order()%r!= 0:
    i+=1
    E_t = EllipticCurve([E.a4(), E.a6() * beta**i])

G2_t = (E_t.order()//r)*E_t.random_point()
while G2_t.is_zero():
    G2_t = (E_t.order()//r)*E_t.random_point()

sqrt3_beta = (X**3-beta**i).roots()[0][0]
sqrt_beta = (beta**i).sqrt()

G2 = E6(G2_t[0] /sqrt3_beta, G2_t[1]/sqrt_beta)
