import numpy as np
import optionpricing

def make_tridiag(n, a, b, c):
    A = np.zeros((n,n))

    A[0,0] = a
    A[0,1] = b

    for i in range(1,n-1,1):
        A[i,i] = a
        A[i,i+1] = b
        A[i,i-1] = c
        
    A[-1,-1] = a
    A[-1,-2] = c
    return A

A = make_tridiag(10,2,1,-1)
x = np.arange(10)

assert np.allclose(np.dot(A, x), optionpricing.tridiag_mult(2.0, 1.0, -1.0, x))
print(np.dot(A, x))
print(optionpricing.tridiag_mult(2.0, 1.0, -1.0, x))

assert np.allclose(np.dot(np.linalg.inv(A), x), optionpricing.tridiag_solve(2.0, 1.0, -1.0, x))
print(np.dot(np.linalg.inv(A), x))
print(optionpricing.tridiag_solve(2.0, 1.0, -1.0, x))
