from sympy import symbols

V_i_j, V_iPlus1_j, V_iMinus1_j, V_i_jPlus1, V_iPlus1_jPlus1, V_iMinus1_jPlus1 = symbols("V_i_j V_iPlus1_j V_iMinus1_j V_i_jPlus1 V_iPlus1_jPlus1 V_iMinus1_jPlus1")
alpha, beta, gamma, dt, dx, c = symbols("alpha beta gamma dt dx c")

lhs = V_i_jPlus1/dt - (c * (alpha * ((V_iPlus1_jPlus1 - 2*V_i_jPlus1 + V_iMinus1_jPlus1)/(dx*dx)) + beta * ((V_iPlus1_jPlus1 - V_iMinus1_jPlus1)/(2*dx)) + gamma * V_i_jPlus1))
rhs = V_i_j/dt + ((1-c) * (alpha * ((V_iPlus1_j - 2*V_i_j + V_iMinus1_j)/(dx*dx)) + beta * ((V_iPlus1_j - V_iMinus1_j)/(2*dx)) + gamma * V_i_j))

print("a1 = ", lhs.subs({V_i_jPlus1:1, V_iPlus1_jPlus1:0, V_iMinus1_jPlus1:0}))
print("b1 = ", lhs.subs({V_i_jPlus1:0, V_iPlus1_jPlus1:1, V_iMinus1_jPlus1:0}))
print("c1 = ", lhs.subs({V_i_jPlus1:0, V_iPlus1_jPlus1:0, V_iMinus1_jPlus1:1}))

print("a0 = ", rhs.subs({V_i_j:1, V_iPlus1_j:0, V_iMinus1_j:0}))
print("b0 = ", rhs.subs({V_i_j:0, V_iPlus1_j:1, V_iMinus1_j:0}))
print("c0 = ", rhs.subs({V_i_j:0, V_iPlus1_j:0, V_iMinus1_j:1}))