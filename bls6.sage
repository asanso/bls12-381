#BLS6 curve taken from https://eprint.iacr.org/2019/431.pdf Table 15

def computeS(n,c):
    """
    (Algorithm 2.35: Lenstra, Verheul: An overview of the XTR public key system)
    Computes S_n(c) as defined in Lenstra, Verheul.
    Parameters:
    (int) n>0;
    (GF(p^2)) c
    Returns:
    (tup) S_n(c)=(c_{n-1},c_n,c_{n+1})
    """
    if type(c) != type(Fp6(1)):
        return(computeS(n,Fp6(c)))
    if n==0:
        return(c**p,Fp6(3),c)
    elif n<0:
        (a1,a2,a3)=computeS(-n,c)
        return(a3**p,a2**p,a1**p)
    elif n==1:
        return(Fp6(3),c,c**2-2*c**p)
    elif n==2:
        (c0,c1,c2)=computeS(1,c)
        return(c1,c2,c*c2-c**p*c1+c0)
    elif n==3:
        (c1,c2,c3)=computeS(2,c)
        return(c2,c3,c*c3-c**p*c2+c1)
    else:
        if n%2==1:
            m=n
        else:
            m=n-1
        k=1
        Sk=computeS(3,c) #Sk=(c_{2k},c_{2k+1},c_{2k+2})

        binario=bin((m-1)/2)[2::]
        r=len(binario)
        contador=0
        if r!= 1:
            for j in range(r):
                if j==0:
                    continue
                if eval(binario[j])==0:
                    Sk=(Sk[0]**2-2*Sk[0]**p,
                         Sk[0]*Sk[1]-c**p*Sk[1]**p+Sk[2]**p,
                         Sk[1]**2-2*Sk[1]**p) #S2k=(c_{4k},c_{4k+1},c_{4k+2})
                    k=2*k
                else:
                    Sk=(Sk[1]**2-2*Sk[1]**p,
                         Sk[2]*Sk[1]-c*Sk[1]**p+Sk[0]**p,
                         Sk[2]**2-2*Sk[2]**p) #S2k+1=(c_{4k+2},c_{4k+3},c_{4k+4})
                    k=2*k+1
        if m==n:
            return(Sk)
        else:
            return(Sk[1],Sk[2],c*Sk[2]-c**p*Sk[1]+Sk[0])


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

a = 3
b = 5

aG1G2 = (a*G1).tate_pairing(G2, r, 6)
bG1G2 = (b*G1).tate_pairing(G2, r, 6)

# standard FF DH
assert aG1G2**b == bG1G2**a

tra = aG1G2 + aG1G2**(p**2) + aG1G2**(p**4) 
trb = bG1G2 + bG1G2**(p**2) + bG1G2**(p**4)  

#XTR DH

assert computeS(a,trb)[1] ==  computeS(b,tra)



