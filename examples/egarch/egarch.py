import numpy as np

def g(eps, theta, lmda):
    return theta * eps + lmda * np.abs(eps)

def log_likelihood(eps, sigma_sq):
    return -np.log(np.pi * sigma_sq) - 0.5 * eps * eps / sigma_sq

def loss(x, deltas):
    print(x)
    alpha, beta, omega, theta, lmda = x

    l = 0.0
    sigma_sq = np.cov(deltas)
    n = len(deltas)

    for eps, epsn in zip(deltas[:-1], deltas[1:]):
        sigma_sq2 = np.exp(omega + beta * g(eps, theta, lmda) + alpha * np.log(sigma_sq))
        l += log_likelihood(epsn, sigma_sq2) / n
        sigma_sq = sigma_sq2

    return -l

def egarch_eval(x, deltas):
    alpha, beta, omega, theta, lmda = x
    
    sigma_sq = np.cov(deltas)
    vol = [sigma_sq]

    for eps in deltas:
        sigma_sq = np.exp(omega + beta * g(eps, theta, lmda) + alpha * np.log(sigma_sq))
        vol.append(sigma_sq)

    return vol