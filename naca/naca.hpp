#include <string>

void naca4_cambered(double m, double p, double t, double c, int n, double xc[],
                    double xu[], double yu[], double xl[], double yl[]);

double *naca4_symmetric(double t, double c, int n, double x[]);
