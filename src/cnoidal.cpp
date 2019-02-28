#include <float.h>
#include <src/cnoidal.h>
#include <stdio.h>
#include <algorithm>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <cmath>

cnoidal::cnoidal() {}

cnoidal::~cnoidal() {}

cnoidal::cnoidal(double H, double T, double d) : H(H), T(T), d(d) {
    calculateL();
    C = sqrt(9.81 * (d + H));
    printf("C = %f\n", C);
    dimensionsless_period = T * sqrt(9.81 / d);
    if (dimensionsless_period < 7) {
        printf("%f<7: Cnoidal might not be valid\n", dimensionsless_period);
    }
    epsilon = H / d;
    ursell = H * pow(L, 2) / pow(d, 3);
    printf("Ursell = %f\n", ursell);
    calculate_m();
    set_cnoidal_velocity();
    set_L_cnoidal();
    printf("m = %f\nc = %f\nL = %f\n", m, c_cnoidal, L_cnoidal);
    alpha = get_alpha(epsilon, m);
}

void cnoidal::calculateL() {
    L = 1.0;
    double wavelength = 10.0;
    while (std::abs(L - wavelength) > 0.0001) {
        wavelength = L;
        L = 9.81 * pow(T, 2) / (2.0 * M_PI) *
            std::tanh(2.0 * M_PI * d / wavelength);
    }

    printf("Wavelength is equal to %f m\n", L);
}

double cnoidal::eta(double x, double t) {
    // third order cnoidal
    double _C = c_cnoidal;
    double eps_over_m = epsilon / m;

    /*return (1.0 + eps_over_m * m * pow(cn(alpha * (x - _C * t) / d, m), 2.0) +
            pow(eps_over_m, 2.0) *
                (-3.0 / 4.0 * pow(m, 2) *
                     pow(cn(alpha * (x - _C * t) / d, m), 2) +
                 3.0 / 4.0 * pow(m, 2.0) *
                     pow(cn(alpha * (x - _C * t) / d, m), 4)) +
            pow(eps_over_m, 3) *
                ((-61.0 / 80.0 * pow(m, 2) + 111.0 / 80.0 * pow(m, 3)) *
                     pow(cn(alpha * (x - _C * t) / d, m), 2) +
                 (61.0 / 80.0 * pow(m, 2) - 53.0 / 20.0 * pow(m, 3)) *
                     pow(cn(alpha * (x - _C * t) / d, m), 4) +
                 101.0 / 80.0 * pow(m, 3) *
                     pow(cn(alpha * (x - _C * t) / d, m), 6))) *
           H;*/
    double Delta = get_L_cnoidal(m) / (2.0 * K(m));

    double eta2 = (H / m) * (1.0 - m - E(m) / K(m));

    return (pow(cn((x - c_cnoidal * t) / Delta, m), 2)) * H + eta2;
}

double cnoidal::get_alpha(double epsilon, double m) {
    return (m < 0.96)
               ? sqrt(3.0 * epsilon / 4.0) *
                     (1.0 - epsilon * (5.0 / 8.0) +
                      (71.0 / 128.0) * pow(epsilon, 2) -
                      (100627.0 / 179200.0) * pow(epsilon, 3) +
                      (16259737.0 / 28672000.0) * pow(epsilon, 4))
               : sqrt(3.0 * epsilon / (4.0 * m)) *
                     (1.0 + epsilon / m * (1.0 / 4.0 - 7.0 / 8.0 * m) +
                      pow(epsilon / m, 2) * (1.0 / 32.0 - 11.0 / 32.0 * m +
                                             111.0 / 128.0 * pow(m, 2)));
}

double cnoidal::K(double m) {
    if (m < 0 || m > 1) printf("Wrong modulus for K\n");
    return boost::math::ellint_1(m);
}
double cnoidal::cn(double z, double m) { return boost::math::jacobi_cn(m, z); }

double cnoidal::K_a(double m) {
    return 2.0 * M_PI / (pow(1.0 + pow(m, 0.25), 2));
}

double cnoidal::q1(double m) { return exp((-1.0) * M_PI * K(m) / K_a(m)); }

double cnoidal::w(double z, double m) { return M_PI * z / (2.0 * K_a(m)); }

double cnoidal::lambda(double m) {
    // third order
    return d * sqrt(16.0 / 3.0 * m * d / H) * K(m);
}

double cnoidal::fac(double n) {
    return (n == 1 || n == 0) ? 1 : n * fac(n - 1);
}

void cnoidal::calculate_m() {
    double m_init = 1.0 - 2 * ds;  // stay away from 1.0
    m = m_init;
    int i = 0;
    double df1, df2, f1, f2, f, df;
    double dx = ds;
    do {
        m = m_init;
        f = T * get_cnoidal_velocity(m) - get_L_cnoidal(m);
        f1 = T * get_cnoidal_velocity(m - dx) - get_L_cnoidal(m - dx);
        f2 = T * get_cnoidal_velocity(m + dx) - get_L_cnoidal(m + dx);
        df = (f2 - f1) / (dx);
        m_init = m - f / df;
        i++;
    } while (abs(f) > ds);
    printf("found m = %f after %d iterations\n", m, i);
}

double cnoidal::E(double m) { return boost::math::ellint_2(m); }

double cnoidal::get_cnoidal_velocity(double m) {
    return c_cnoidal =
               sqrt(d * 9.81) *
               (1.0 + (H / (m * d)) * (1.0 - 0.5 * m - 1.5 * (E(m) / K(m))));
}

void cnoidal::set_cnoidal_velocity() { c_cnoidal = get_cnoidal_velocity(m); }

void cnoidal::set_L_cnoidal() { L_cnoidal = get_L_cnoidal(m); }

double cnoidal::get_L_cnoidal(double m) {
    return d *
           sqrt((16.0 / 3.0) * (m * d / H) * get_cnoidal_velocity(m) /
                sqrt(9.81 * d)) *
           K(m);
}
