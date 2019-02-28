#include <src/stokesIV.h>
#include <stdio.h>
#include <cmath>
StokesIV::StokesIV() {}

StokesIV::~StokesIV() {}

StokesIV::StokesIV(double H, double T, double d) : H(H), T(T), d(d) {
    calculateL();
    omega = 2 * M_PI / T;
}

void StokesIV::calculateL() {
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

double StokesIV::eta(double x, double t) {
    double theta = k * x - omega * t;
    double a = H / 2.0;
    double third_order =
        a * ((1.0 - 1.0 / 16.0 * (k * a) * (k * a)) * cos(theta) +
             0.5 * (k * a) * cos(2 * theta) +
             3.0 / 8.0 * pow(k * a, 2) * cos(3 * theta));

    return third_order;
}
