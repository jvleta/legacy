#include <string>
#include <vector>

void naca4_cambered(double m, double p, double t, double c, int n, const std::vector<double>& xc,
                    std::vector<double>& xu, std::vector<double>& yu, std::vector<double>& xl, std::vector<double>& yl);

std::vector<double> naca4_symmetric(double t, double c, int n, const std::vector<double>& x);
