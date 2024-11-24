import numpy as np
import matplotlib.pyplot as plt

import optionpricing
from tridiag import make_tridiag
from black_scholes import call

alpha = 0.025
beta = 0.0
gamma = 0.0
pde = optionpricing.PDESolver(alpha, beta, gamma)

dx = 0.02
x = np.arange(-4.0, 4.0+dx, dx)
pdf = x.copy()
pdf[pdf<0.0] = 0.0

y = pdf.copy()

dt = 0.1
c = 0.5

for i in range(500):
    y = pde.step(y, c, dt, dx)

plt.plot(x, pdf)
plt.plot(x, y)
plt.savefig("examples/pde.png")