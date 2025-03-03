#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "naca.hpp"

void naca4_cambered(double m, double p, double t, double c, int n, const std::vector<double>& xc,
                    std::vector<double>& xu, std::vector<double>& yu, std::vector<double>& xl, std::vector<double>& yl) {
  //
  //  Purpose:
  //
  //    NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
  //
  //  Licensing:
  //
  //    This code is distributed under the MIT license.
  //
  //  Modified:
  //
  //    22 May 2014
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
  //    "The characteristics of 78 related airfoil sections from tests in
  //    the variable-density wind tunnel",
  //    NACA Report 460, 1933.
  //
  //  Parameters:
  //
  //    Input, double M, the maximum camber.
  //    0.0 < M.
  //
  //    Input, double P, the location of maximum camber.
  //    0.0 < P < 1.0
  //
  //    Input, double T, the maximum relative thickness.
  //    0.0 < T <= 1.0
  //
  //    Input, double C, the chord length.
  //    0.0 < C.
  //
  //    Input, int N, the number of sample points.
  //
  //    Input, double XC[N], points along the chord length.
  //    0.0 <= XC(*) <= C.
  //
  //    Output, double XU[N], YU[N], XL[N], YL[N], for each value of
  //    XC, measured along the camber line, the corresponding values (XU,YU)
  //    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
  //

  for (int i = 0; i < n; i++) {

    double divisor = [&]() {
      if (0.0 <= xc[i] / c && xc[i] / c <= p) {
        return p * p;
      } else if (p <= xc[i] / c && xc[i] / c <= 1.0) {
        return pow(1.0 - p, 2);
      } else {
        return 1.0;
      }
    }();

    double dycdx = 2.0 * m * (p - xc[i] / c) / divisor;
    double theta = std::atan(dycdx);

    double yt = 5.0 * t * c *
                (0.2969 * std::sqrt(xc[i] / c) +
                 ((((-0.1015) * (xc[i] / c) + 0.2843) * (xc[i] / c) - 0.3516) *
                      (xc[i] / c) -
                  0.1260) *
                     (xc[i] / c));

    double yc = [&]() {
      if (0.0 <= xc[i] / c && xc[i] / c <= p) {
        return m * xc[i] * (2.0 * p - xc[i] / c) / p / p;
      } else if (p <= xc[i] / c && xc[i] / c <= 1.0) {
        return m * (xc[i] - c) * (2.0 * p - xc[i] / c - 1.0) / (1.0 - p) /
               (1.0 - p);
      } else {
        return 0.0;
      }
    }();

    xu[i] = xc[i] - yt * std::sin(theta);
    yu[i] = yc + yt * std::cos(theta);
    xl[i] = xc[i] + yt * std::sin(theta);
    yl[i] = yc - yt * std::cos(theta);
  }
  return;
}

std::vector<double> naca4_symmetric(double t, double c, int n, const std::vector<double>& x) {
  //
  //  Purpose:
  //
  //    NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
  //
  //  Licensing:
  //
  //    This code is distributed under the MIT license.
  //
  //  Modified:
  //
  //    22 May 2014
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
  //    "The characteristics of 78 related airfoil sections from tests in
  //    the variable-density wind tunnel",
  //    NACA Report 460, 1933.
  //
  //  Parameters:
  //
  //    Input, double T, the maximum relative thickness.
  //
  //    Input, double C, the chord length.
  //
  //    Input, int N, the number of sample points.
  //
  //    Input, double X[N], points along the chord length.
  //    0.0 <= X(*) <= C.
  //
  //    Output, double NACA4_SYMMETRIC[N], for each value of X, the
  //    corresponding value of Y so that (X,Y) is on the upper wing surface, and
  //    (X,-Y) is on the lower wing surface.
  //
  std::vector<double> y(n);

  for (int i = 0; i < n; i++) {
    y[i] = 5.0 * t * c *
           (0.2969 * std::sqrt(x[i] / c) +
            ((((-0.1015) * (x[i] / c) + 0.2843) * (x[i] / c) - 0.3516) *
                 (x[i] / c) -
             0.1260) *
                (x[i] / c));
  }

  return y;
}
