#include <src/stokesIV.h>
#include <stdio.h>
#include <cmath>
StokesIV::StokesIV() {}

StokesIV::~StokesIV() {}

StokesIV::StokesIV(double H, double T, double d) : H(H), T(T), d(d) {
    calculateL();
    omega = 2 * M_PI / T;
    calculate_coefficients();
    epsilon = k * H / 2.0;
}

void StokesIV::calculateL() {
    L = 1;
    double wavelength = 10;
    while (std::abs(L - wavelength) > 0.0001) {
        wavelength = L;
        L = 9.81 * pow(T, 2) / (2.0 * M_PI) *
            std::tanh(2.0 * M_PI * d / wavelength);
    }
    k = 2.0 * M_PI / L;
    C = L / T;
    printf("Wavelength is equal to %f m\n", L);
}

double StokesIV::eta(double x, double t) {
    double theta = k * x - omega * t;
    double a = H / 2.0;
    double third_order =
        a * ((1.0 - 1.0 / 16.0 * (k * a) * (k * a)) * cos(theta) +
             0.5 * (k * a) * cos(2 * theta) +
             3.0 / 8.0 * pow(k * a, 2) * cos(3 * theta));

    // double etaFifth[5];
    // for (int i = 0; i < order; i++) {
    //    etaFifth[i] = 0;
    //    for (int j = 0; j <= i; j++) {
    //        etaFifth[i] += B[i][j] * cos((double)(j + 1) * k * (x - C * t)) *
    //                       pow(epsilon, i + 1);
    //    }
    //}
    // double eta_fi = 0;
    // for (int i = 0; i < order; i++) {
    //    eta_fi += etaFifth[i];
    //}
    // eta_fi = 0;
    // double c = 1.0;
    // for (int i = 0; i < N; i++) {
    //    if (i == 0 || (i == (N - 1))) c = 0.5;
    //    eta_fi += Y((double)(i + 1)) * cos((i + 1.) * k * (x - C * t)) * c;
    //    c = 1.0;
    //}
    return third_order;
}

void StokesIV::calculate_coefficients() {
    k = 2 * M_PI / L;
    double kd = k * d;
    double S = 1.0 / cosh(2 * kd);
    B[0][0] = 1.0;
    B[1][1] = 1.0 / (tanh(kd)) * (1 + 2. * S) / (2. * (1. - S));
    B[2][0] = -3.0 * (1.0 + 3. * S + 3. * S * S + 2. * S * S * S) /
              (8. * pow(1. - S, 3));
    B[2][2] = (-1.0) * B[2][0];
    B[3][1] = 1.0 / tanh(kd) *
              (6.0 - 26.0 * S - 182.0 * pow(S, 2) - 204.0 * pow(S, 3) -
               25.0 * pow(S, 4) + 26.0 * pow(S, 5)) /
              (6.0 * (3.0 + 2.0 * S) * pow(1 - S, 4));
    B[3][3] = 1.0 / tanh(kd) *
              (24. + 92. * S + 122. * pow(S, 2) + 66. * pow(S, 3) +
               67. * pow(S, 4) + 34. * pow(S, 5)) /
              (24. * (3.0 + 2.0 * S) * pow(1 - S, 4));
    B[4][0] = 9 *
              (132. + 17 * S - 2216. * pow(S, 2) - 5897. * pow(S, 3) -
               6292. * pow(S, 4) - 2687. * pow(S, 5) + 194. * pow(S, 6) +
               467. * pow(S, 7) + 82. * pow(S, 8)) /
              (128. * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
    B[4][4] = 5.0 *
              (300. + 1579. * C + 3176. * pow(C, 2) + 2949. * pow(C, 3) +
               1188. * pow(C, 4) + 675. * pow(C, 5) + 1326. * pow(C, 6) +
               827. * pow(C, 7) + 130. * pow(C, 8)) /
              (384. * (3. + 2. * S) * (4. + S) * pow(1 - S, 6));
    B[4][0] = (-1.0) * (B[4][2] + B[4][4]);
}

double StokesIV::Y(double j) {
    double I = 0;
    double eta_m = 0;
    double c = 1.0;
    for (int i = 0; i < N; i++) {
        eta_m = d + 0.5 * H * cos((double)(i)*M_PI / order);
        if (i == 0 || (i == (N - 1))) c = 0.5;
        I += (2.0 / order) * k * eta_m *
             cos((double)(i)*j * M_PI / ((double)(N))) * c;
        c = 1.0;
    }
    return I;
}
