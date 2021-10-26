q, p = 2**23, 2**19
k, l = 3, 2
n = 256
h = 60
eta = q // (2 * p)
gamma1 = q // 16


# def pm_mod(a, q):
#     result = QQ(a) % q
#     if result > q // 2:
#         result -= q
#     return result

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
    return vector([ring([(ring_mod/prev_mod * QQ(coeff)).round('up')
                         for coeff in poly])
                   for poly in vec])


def simple_round_vec(ring, vec):
    return vector([ring([coeff.round('up') for coeff in poly])
                   for poly in vec])


def coerce_val(ring, prev_mod, poly):
    return ring([coeff for coeff in poly])


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
P_cont.<X> = QQ[]
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
xi2 = (simple_round_vec(R_cont, coerce_val(R_cont, p, c) * s2) -
       coerce_val(R_cont, p, c) * s2)
nu = simple_round_vec(R, xi1 - xi2)

lhs = round_vec(Rp, p, q, A * z) - coerce_val(Rp, p, c) * t
rhs = (w - simple_round_vec(Rp, coerce_val(R_cont, p, c) * s2) -
       coerce_vec(Rp, p, nu))


# helper function for debugging
def err_pq(vec):
    return round_vec(R_cont, p, q, vec) - p/q * coerce_vec(R_cont, q, vec)

err_pqAy = err_pq(A * y)
err_cpqAs1 = err_pq(A * coerce_val(Rq, q, c) * s1)

lhs1 = (round_vec(Rp, p, q, A * y) +
        round_vec(Rp, p, q, A * coerce_val(Rq, q, c) * s1) +
        simple_round_vec(Rp, -err_pqAy - err_cpqAs1) -
        coerce_val(Rp, p, c) * t)

print(lhs - lhs1)

'''
# okay... make it simpler
should_be_zero = (round_vec(Rp, p, q, A * z) -
                  (round_vec(Rp, p, q, A * y) +
                   round_vec(Rp, p, q, A * coerce_val(Rq, q, c) * s1) +
                   simple_round_vec(Rp, -err_pqAy - err_cpqAs1)))


print(list(p/q * coerce_vec(R_cont, q, (A * z))[0])[0])
print(list(p/q * coerce_vec(R_cont, q, (A * y))[0])[0])
print(list(p/q * coerce_vec(
    R_cont, q, (A * coerce_val(Rq, q, c) * s1))[0])[0])

print('--')
print(list(round_vec(Rp, p, q, A * z)[0])[0])
print(list(round_vec(Rp, p, q, A * y)[0])[0])
print(list(round_vec(Rp, p, q, A * coerce_val(Rq, q, c) * s1)[0])[0])
print(list(simple_round_vec(Rp, -err_pqAy - err_cpqAs1)[0])[0])

print('--')
print(list((-err_pqAy - err_cpqAs1)[0])[0])
print((list((-err_pqAy - err_cpqAs1)[0])[0]).round('up'))

print('--')
print(list(should_be_zero[0]))
'''
