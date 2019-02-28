#pragma once
class StokesII {
   public:
    double H, T, d;
    double L;
    double k = 0;
    double omega = 0;
    StokesII();
    StokesII(double H, double T, double d);
    ~StokesII();

    double eta(double x, double t);
    void calculateL();
};
