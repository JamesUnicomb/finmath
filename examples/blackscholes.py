import numpy as np
from math import erf

def N(x):
    return 1.0/2.0 + 1.0 / 2.0 * erf(x / np.sqrt(2.0))

def d1(S, K, r, q, sigma, T):
    return 1.0 / (sigma * np.sqrt(T)) * (np.log(S / K) + (r - q + 1.0 / 2.0 * sigma * sigma) * (T))

def d2(S, K, r, q, sigma, T):
    return 1.0 / (sigma * np.sqrt(T)) * (np.log(S / K) + (r - q - 1.0 / 2.0 * sigma * sigma) * (T))

def C(S, K, r, q, sigma, T):
    return S * np.exp(-q * T) * N(d1(S, K, r, q, sigma, T)) - K * np.exp(-r * T) * N(d2(S, K, r, q, sigma, T))

def P(S, K, r, q, sigma, T):
    return K * np.exp(-r * T) * N(-d2(S, K, r, q, sigma, T)) - S * np.exp(-q * T) * N(-d1(S, K, r, q, sigma, T))

def l(r, q, sigma):
    return (r - q + np.square(sigma)/2.0) / np.square(sigma)

def y(r, q, H, S, K, sigma, T):
    return np.log(np.square(H) / (S*K)) / (sigma * np.sqrt(T)) + l(r, q, sigma) * sigma * np.sqrt(T)

def x1(r, q, H, S, K, sigma, T):
    return np.log(S / H) / (sigma * np.sqrt(T)) + l(r, q, sigma) * sigma * np.sqrt(T)

def y1(r, q, H, S, K, sigma, T):
    return np.log(H / S) / (sigma * np.sqrt(T)) + l(r, q, sigma) * sigma * np.sqrt(T)

def Cui(S, K, H, r, q, sigma, T):
    return S * np.exp(-q * T) * N(x1(r, q, H, S, K, sigma, T)) - \
        K * np.exp(-r * T) * N(x1(r, q, H, S, K, sigma, T) - sigma * np.sqrt(T)) - \
            S * np.exp(-q * T) * np.power(H / S, 2.0 * l(r, q, sigma)) * (N(-y(r, q, H, S, K, sigma, T)) - N(-y1(r, q, H, S, K, sigma, T))) + \
                K * np.exp(-r * T) * np.power(H / S, 2.0 * l(r, q, sigma) - 2.0) * (N(-y(r, q, H, S, K, sigma, T)+ sigma * np.sqrt(T)) - N(-y1(r, q, H, S, K, sigma, T) + sigma * np.sqrt(T)))

def Cuo(S, K, H, r, q, sigma, T):
    return C(S, K, r, q, sigma, T) - Cui(S, K, H, r, q, sigma, T)