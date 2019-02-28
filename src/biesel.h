#pragma once
#include <vector>

class Biesel {
   public:
    double H, h, T;
    double L, m_k;
    double S0;
    int N = 19;
    double ds = 0.0001;

    std::vector<double> k;
    Biesel(double h, double H, double T);
    ~Biesel();

    void calculateL();

    void calculateS0();

    void calculatek(int n);

    double e_piston_z(double z);

    double omega(double n);

    double e_t(double e_z, double t, double n);

    double c_piston(int n);

    double integrate_ezcosh(double under, double upper, int n);

    double eta(double x, double t);

    double k_function(double x);

    double k_function_der(double x);
};
