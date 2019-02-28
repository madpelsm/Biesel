#include <src/stokesII.h>
#include <stdio.h>
#include <cmath>
StokesII::StokesII() {}

StokesII::~StokesII() {}

StokesII::StokesII(double H, double T, double d) : H(H), T(T), d(d) {
    calculateL();
    omega = 2 * M_PI / T;
}

void StokesII::calculateL() {
    L = 1;
    double wavelength = 10;
    while (std::abs(L - wavelength) > 0.0001) {
        wavelength = L;
        L = 9.81 * pow(T, 2) / (2 * M_PI) *
            std::tanh(2 * M_PI * d / wavelength);
    }
    k = 2 * M_PI / L;

    printf("Wavelength is equal to %f m\n", L);
}

double StokesII::eta(double x, double t) {
    double a1 = H / 2;
    double epsilon = 2 * M_PI / L * a1;
    double theta = k * x - omega * t;

    double a2 =
        a1 / 4 * cosh(k * d) / pow(sinh(k * d), 3) * (2 + cosh(2 * k * d));
    return a1 * cos(theta) + epsilon * a2 * cos(2 * theta);
}
