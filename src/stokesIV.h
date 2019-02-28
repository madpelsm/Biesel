class StokesIV {
   public:
    double d, H, T, L;
    double k = 0;
    double omega = 0;
    StokesIV();
    StokesIV(double H, double T, double d);
    ~StokesIV();

    double eta(double x, double t);
    void calculateL();

};
