#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <vector>

void tridiag_mult(double a, double b, double c, std::vector<double> &r, std::vector<double> &u)
{
    int j;    
    int n = r.size();

    u[0] = a*r[0] + b*r[1];

    for (j = 1; j < n-1; j++)
        u[j] = c*r[j-1] + a*r[j] + b*r[j+1];

    u[n-1] = c*r[n-2] + a*r[n-1];
}

std::vector<double> tridiag_mult(double a, double b, double c, std::vector<double> &r)
{
    int n = r.size();

    std::vector<double> u(n, 0.0);

    tridiag_mult(a, b, c, r, u);

    return u;
}

void tridiag_solve(double a, double b, double c, std::vector<double> &r, std::vector<double> &u)
{
    int j;
    double bet;
    
    int n = r.size();
    std::vector<double> gam(n);

    u[0] = r[0] / (bet = a);

    for (j = 1; j < n; j++)
    {
        gam[j] = b / bet;
        bet = a - c * gam[j];
        u[j] = (r[j] - c * u[j - 1]) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
        u[j] -= gam[j + 1] * u[j + 1];
}

std::vector<double> tridiag_solve(double a, double b, double c, std::vector<double> &r)
{
    int n = r.size();

    std::vector<double> u(n, 0.0);

    tridiag_solve(a, b, c, r, u);

    return u;
}

#endif