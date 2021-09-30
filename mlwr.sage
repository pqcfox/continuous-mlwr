q, p = 2**23, 2**19
k, l = 3, 2
n = 256
h = 60
eta = q // (2 * p)
gamma1 = q // 16


def pm_mod(a, q):
    result = a % q
    if result > q // 2:
        result -= q
    return result


def random_mat(ring, ring_mod, ring_deg, nrows, ncols):
    return matrix([[ring([randrange(0, ring_mod)
                          for _ in range(ring_deg)])
                    for _ in range(ncols)]
                   for _ in range(nrows)])


def uniform_vec(ring, val_range, size):
    return vector([ring([randint(-val_range, val_range)
                      for _ in range(n)])
                   for _ in range(size)])


def round_vec(ring, ring_mod, prev_mod, vec):
    return vector([ring([round(ring_mod/prev_mod *
                               pm_mod(int(coeff), prev_mod))
                         for coeff in poly])
                   for poly in vec])


def coerce_val(ring, prev_mod, poly):
    return ring([pm_mod(int(coeff), prev_mod) for coeff in poly])


def coerce_vec(ring, prev_mod, vec):
    return vector([coerce_val(ring, prev_mod, poly)
                   for poly in vec])


def gen_chal(ring, ring_deg, num_hot):
    coeffs = [0 for _ in range(ring_deg)]
    for index in sample(range(len(coeffs)), num_hot):
        coeffs[index] = choice([-1, 1])
    return ring(coeffs)


P.<X> = ZZ[]
Pq.<X> = Integers(q)[]
Pp.<X> = Integers(p)[]
P_cont.<X> = RR[]
R.<x> = QuotientRing(P, X^n + 1)
Rq.<x> = QuotientRing(Pq, X^n + 1)
Rp.<x> = QuotientRing(Pp, X^n + 1)
R_cont.<x> = QuotientRing(P_cont, X^n + 1)

A = random_mat(Rq, q, n, k, l)
s1 = uniform_vec(Rq, eta, l)
t = round_vec(Rp, p, q, A * s1)

y = uniform_vec(Rq, gamma1 - 1, l)
w = round_vec(Rp, p, q, A * y)

c = gen_chal(R, n, h)
s2 = coerce_vec(R_cont, p, t) - p/q * coerce_vec(R_cont, q, A * s1)
z = y + coerce_val(Rq, q, c) * s1
xi1 = round_vec(R_cont, p, q, A * y) - p/q * coerce_vec(R_cont, q, A * y)
xi2 = (round_vec(R_cont, p, p, coerce_val(R_cont, p, c) * s2) -
       coerce_val(R_cont, p, c) * s2)
nu = round_vec(R, p, p, xi1 - xi2)

lhs = round_vec(Rp, p, q, A * z) - coerce_val(Rp, p, c) * t
rhs = (w - round_vec(Rp, p, p, coerce_val(R_cont, p, c) * s2) -
       coerce_vec(Rp, p, nu))

print(lhs - rhs)  # should be zero
