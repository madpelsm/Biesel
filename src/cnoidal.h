#pragma once

class cnoidal {
   public:
    double H, T, d, C, L;
    double trough_elevation;
    double m = 0;
    double dimensionsless_period;
    double c_cnoidal;
    double L_cnoidal;
    double epsilon;
    double alpha;
    double ursell;
    int N = 60;
    double ds = 1.0e-9;

    cnoidal();
    ~cnoidal();
    cnoidal(double H, double T, double d);
    double eta(double x, double t);

    void calculateL();
    double K(double m);
    double cn(double z, double m);
    double K_a(double m);
    double q1(double m);
    double w(double z, double m);
    double lambda(double m);
    void calculate_m();
    double E(double m);
    void set_cnoidal_velocity();
    void set_L_cnoidal();
    double get_alpha(double epsilon, double m);
    double get_L_cnoidal(double m);
    double fac(double n);
    double get_cnoidal_velocity(double m);
};
