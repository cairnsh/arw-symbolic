from sympy import *
LAMBDA = lambda m: symbols("λ_%d" % m)
FRAC = lambda p, q: Rational(p, q)
STOP = -999

def polynomialize(z):
    z = [simplify(j) for j in z]
    the_lcm = lcm_list([denom(j) for j in z])
    z = [cancel(j * the_lcm) for j in z]
    the_gcd = gcd_list(z)
    z = [cancel(j / the_gcd) for j in z]
    rr = Matrix(z)
    return rr