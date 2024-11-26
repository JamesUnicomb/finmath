#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <vector>

void tridiag_solve(double a, double b, double c, std::vector<double> &r, std::vector<double> &u) {
    int j;
	double bet;

	int n=u.size();
	std::vector<double> gam(n,0.0);
	if (b == 0.0) {
        throw std::invalid_argument("Error 1 in tridag");
    }
    bet=b;
	u[0]=r[0]/bet;
	for (j=1;j<n;j++) {
		gam[j]=c/bet;
		bet=b-a*gam[j];
		if (bet == 0.0) {
            throw std::invalid_argument("Error 2 in tridag");
        }
		u[j]=(r[j]-a*u[j-1])/bet;
	}
	for (j=(n-2);j>=0;j--)
		u[j] -= gam[j+1]*u[j+1];
}

std::vector<double> tridiag_solve(double a, double b, double c, std::vector<double> &r) {
    std::vector<double> u(r.size(), 0.0);
    tridiag_solve(a, b, c, r, u);
    return u;
}

void tridiag_mult(double a, double b, double c, std::vector<double> &r, std::vector<double> &u) {
    int j;
	double bet;

	int n=u.size();
	
    u[0] = b*r[0] + c*r[1];
    u[n-1] = a*r[n-2] + b*r[n-1];
    for (j = 1; j < n-1; j++) {
        u[j] = a*r[j-1] + b*r[j] + c*r[j+1];
    }
}

std::vector<double> tridiag_mult(double a, double b, double c, std::vector<double> &r) {
    std::vector<double> u(r.size(), 0.0);
    tridiag_mult(a, b, c, r, u);
    return u;
}

#endif