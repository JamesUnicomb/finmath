import numpy as np

def ewma_eval(alpha, deltas):
    vol = [0.0]

    for eps in deltas:
        vol.append((1.0 - alpha) * eps * eps + alpha * vol[-1])

    return vol