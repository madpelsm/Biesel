class StokesIV {
   public:
    double d, H, T, L;
    double k = 0;
    double omega = 0;
    double C = 0;
    double epsilon = 0;
    double B[5][5];
    int order = 5;
    int N = 10;
    StokesIV();
    StokesIV(double H, double T, double d);
    ~StokesIV();

    double eta(double x, double t);

    void calculateL();

    void calculate_coefficients();

    double Y(double j);
};
