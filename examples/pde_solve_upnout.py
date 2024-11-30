import numpy as np
import matplotlib.pyplot as plt

from blackscholes import Cuo
import optionpricing

K = 100
H = 200
r = 0.05
q = 0.0
sigma = 0.24
T = 1.0
L = np.log(H / K)

# du / dt = alpha * beta^2u / dx^2 + beta * du / dx + gamma * u
gamma = -r
beta = (r - q) - 1.0/2.0 * np.square(sigma)
alpha = 1.0 / 2.0 * np.square(sigma)

nx = 100
dt = 0.05
# end points are defined so start index at 1
x = np.linspace(-L, L, nx)[1:-1]
dx = 2*L/nx
print(np.mean(np.diff(x)), dx)
t = np.arange(0.0, T, dt)

print("dt / (dx)^2 = ", dt / (dx * dx))

u = K*np.max(np.column_stack([np.zeros_like(x), np.exp(x)-1.0]), axis=1)

plt.plot(K*np.exp(x), u)

a1 = -0.5 * (alpha * dt / (dx * dx) - beta * dt / (2.0 * dx))
b1 = 1.0 - gamma * dt / 2.0 + alpha * dt / (dx * dx)
c1 = -0.5 * (alpha * dt / (dx * dx) + beta * dt / (2.0 * dx))

a2 = 0.5 * (alpha * dt / (dx * dx) - beta * dt / (2.0 * dx))
b2 = 1.0 + gamma * dt / 2.0 - alpha * dt / (dx * dx)
c2 = 0.5 * (alpha * dt / (dx * dx) + beta * dt / (2.0 * dx))

def phi(tau):
    return 0.0

def psi(tau):
    return 0.0

def boundary(tau):
    r = np.zeros_like(x)
    r[0] = -a1 * phi(tau + dt) + a2 * phi(tau)
    r[-1] = -c1 * psi(tau + dt) + c2 * psi(tau)
    return r

for tau in t:
    u = optionpricing.tridiag_mult(a2, b2, c2, u) + boundary(tau)
    u = optionpricing.tridiag_solve(a1, b1, c1, u)

plt.plot(K*np.exp(x), K*np.max(np.column_stack([np.zeros_like(x), np.exp(x)-1.0]), axis=1), '--')
plt.plot(K*np.exp(x), u)
plt.plot(K*np.exp(x), [Cuo(s, K, H, r, q, sigma, T) for s in K*np.exp(x)], '--')
plt.show()
