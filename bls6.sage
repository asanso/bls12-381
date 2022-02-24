#BLS6 curve taken from https://eprint.iacr.org/2019/431.pdf Table 15

x = -0xdffffffffffffffffbffffffffffffffffffffdfffffffc00001 
r = x^2 - x +1   
p = 1/3*(x-1)^2*(x^2 -x+1) +x 
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
D = fundamental_discriminant(-3)
H = FpX(hilbert_class_polynomial(D))
(j,_) = H.roots()
