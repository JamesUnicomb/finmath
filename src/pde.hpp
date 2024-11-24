#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include "tridiag.hpp"

class PDESolver {
public:
    PDESolver(double alpha, double beta, double gamma) : alpha(alpha), beta(beta), gamma(gamma) {}
    
    std::vector<double> step(std::vector<double> v, double c, double dt, double dx, double r, double S, double K, double currT) {
        // set up crank nicholson constants
        a1 =  -c*(-2*alpha/(dx*dx) + gamma) + 1/dt;
        b1 =  -c*(alpha/(dx*dx) + beta/(2*dx));
        c1 =  -c*(alpha/(dx*dx) - beta/(2*dx));
        a0 =  (1.0-c)*(-2*alpha/(dx*dx) + gamma) + 1/dt;
        b0 =  (1.0-c)*(alpha/(dx*dx) + beta/(2*dx));
        c0 =  (1.0-c)*(alpha/(dx*dx) - beta/(2*dx));

        int n = v.size();
        std::vector<double> tmp(n, 0.0);

        tridiag_mult(a0, b0, c0, v, tmp);
        tmp[n-1] += c0*(S - K*exp(-r*currT)) - c1*(S - K*exp(-r*(currT+dt)));
        tridiag_solve(a1, b1, c1, tmp, v);

        return v;
    }
    double alpha, beta, gamma;
    double a1, b1, c1, a0, b0, c0;
};

#endif