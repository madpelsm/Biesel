#include <src/biesel.h>
#include <stdio.h>
#include <cmath>

Biesel::Biesel(double _h, double _H, double _T) : H(_H), h(_h), T(_T) {
    calculateL();
    calculateS0();
    k.reserve(N);
    for (int n = 0; n < N; n++) {
        calculatek(n);
        calc_piston(n);
    }
}

Biesel::~Biesel() {}

void Biesel::calculateL() {
    L = 1;
    double wavelength = 10;
    while (std::abs(L - wavelength) > 0.0001) {
        wavelength = L;
        L = 9.81 * pow(T, 2) / (2 * M_PI) *
            std::tanh(2 * M_PI * h / wavelength);
    }
    m_k = (2.0 * M_PI) / L;
    printf("Wavelength is equal to %f m\n", L);
}

void Biesel::calculateS0() {
    S0 = (H / (2.0 * pow(sinh(m_k * h), 2))) *
         (sinh(m_k * h) * cosh(m_k * h) + m_k * h);
    printf("S0 is equal to: %fm\n", S0);
}

void Biesel::calculatek(int n) {
    // dirt way of solving it
    if (n == 0)
        k.push_back(m_k);
    else {
        double kn, knn, f, df;
        knn = (double)(n) * (M_PI / h);
        do {
            kn = knn;
            f = k_function(kn);
            df = k_function_der(kn);
            knn = kn - f / df;
        } while (abs(f) >= ds);
        k.push_back(kn);
    }
}

double Biesel::k_function(double x) {
    return x * (-9.81) * tan(x * h) - pow(omega(0), 2);
}

double Biesel::k_function_der(double x) {
    return (-9.81) * tan(x * h) - 9.81 * x * h * (1 + pow(tan(x * h), 2));
}

double Biesel::e_piston_z(double z) { return S0 * 0.5; }

double Biesel::omega(double n) {
    if (n == 0) {
        return 2 * M_PI / T;
    } else {
        return 2 * M_PI * n / T;
    }
}

double Biesel::e_t(double e_z, double t, double n = 1) {
    return e_z * sin(omega(0) * t);
}

void Biesel::calc_piston(int n) {
    double noemer = (sinh(k[n] * h) * cosh(k[n] * h) + k[n] * h);
    double teller = (2.0 * k[n] * integrate_ezcosh(0, h, n));
    double ccc = teller / noemer;
    if (ccc / H < ds) {
        N = n;
        printf("Only need %d terms\n", N);
    }
    c_pist.push_back(ccc);
}

double Biesel::integrate_ezcosh(double under, double upper, int n) {
    double dz = ds;
    double I = 0;
    double z = under;

    while (z < upper) {
        I += (e_piston_z(z - h) * cosh(k[n] * (z))) * dz;
        z += dz;
    }

    //    printf("Integration yielded %f\n", I);

    return I;
}

double Biesel::eta(double x, double t) {
    double ss = H / 2 * cos(omega(0) * t - m_k * x);
    double T = 0;
    for (int n = 1; n < N; n++) {
        T = exp(((-1.0) * k[n]) * x) *
            (sin(omega(0) * t) * (c_pist[n] * sin(k[n] * h)));
        ss += T;
    }
    return ss;
}
