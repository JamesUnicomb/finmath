import numpy as np
import matplotlib.pyplot as plt

import optionpricing
from black_scholes import call

sigma = 0.5
r = 0.2
K = 100
T = 0.25

alpha = 0.5*np.square(sigma)
beta = 0.0
gamma = r
pde = optionpricing.PDESolver(alpha, beta, gamma)

dS = 0.01
S = np.exp(np.arange(np.log(40), np.log(160), dS))
V = S - K
V[V<0.0] = 0.0

v = V.copy()

N = 100
dt = T / N
c = 1.0
t = 0.0

for i in range(N):
    v = pde.step(v, c, dt, dS, r, S[-1], K, t)
    t += dt

print(S)
p = np.array([call(s, K, T, r, sigma) for s in S])
plt.plot(S, V)
plt.plot(S, v)
plt.plot(S, p)
plt.savefig("examples/call.png")