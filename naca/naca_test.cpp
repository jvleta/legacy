#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#include "naca.hpp"

namespace {
  void r8mat_write(std::string output_filename, int m, int n, const std::vector<double>& table) {
    std::ofstream output;
    output.open(output_filename.c_str());

    if (!output) {
      std::cerr << "\n";
      std::cerr << "R8MAT_WRITE - Fatal error!\n";
      std::cerr << "  Could not open the output file.\n";
      std::exit(1);
    }

    for (auto j = 0; j < n; j++) {
      for (auto i = 0; i < m; i++) {
        output << "  " << std::setw(24) << std::setprecision(16)
               << table[i + j * m];
      }
      output << "\n";
    }

    output.close();

    return;
  }

  std::vector<double> r8vec_linspace_new(int n, double a_first, double a_last) {
    std::vector<double> a(n);

    if (n == 1) {
      a[0] = (a_first + a_last) / 2.0;
    } else {
      for (int i = 0; i < n; i++) {
        a[i] = (static_cast<double>(n - 1 - i) * a_first + static_cast<double>(i) * a_last) /
               static_cast<double>(n - 1);
      }
    }
    return a;
  }

  double r8vec_max(int n, const std::vector<double>& r8vec) {
    double value = r8vec[0];

    for (int i = 1; i < n; i++) {
      if (value < r8vec[i]) {
        value = r8vec[i];
      }
    }
    return value;
  }

  double r8vec_min(int n, const std::vector<double>& r8vec) {
    double value = r8vec[0];

    for (int i = 1; i < n; i++) {
      if (r8vec[i] < value) {
        value = r8vec[i];
      }
    }
    return value;
  }
}

void test01() {
  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TEST01 tests NACA4_SYMMETRIC.
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

  std::string command_filename = "symmetric_commands.txt";
  std::ofstream command_unit;
  std::string data_filename = "symmetric_data.txt";
  int n = 51;
  std::vector<double> xy;
  
  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  NACA4_SYMMETRIC evaluates y(x) for a NACA\n";
  std::cout << "  symmetric airfoil defined by a 4-digit code.\n";

  const auto c = 10.0;
  const auto t = 0.15;
  const auto x = r8vec_linspace_new(n, 0.0, c);
  const auto y = naca4_symmetric(t, c, n, x);
  //
  //  Reorganize data into a single object.
  //
  xy.resize(2 * 2 * n);

  for (auto i = 0; i < n; i++) {
    xy[0 + i * 2] = x[i];
    xy[1 + i * 2] = -y[i];
  }
  for (auto i = 0; i < n; i++) {
    xy[0 + (n + i) * 2] = x[n - 1 - i];
    xy[1 + (n + i) * 2] = y[n - 1 - i];
  }
  //
  //  Determine size ratio.
  //
  const auto x_min = r8vec_min(n, x);
  const auto x_max = r8vec_max(n, x);
  const auto y_max = r8vec_max(n, y);
  const auto y_min = -y_max;
  const auto ratio = (y_max - y_min) / (x_max - x_min);
  //
  //  Save data to a file.
  //
  r8mat_write(data_filename, 2, 2 * n, xy);
  std::cout << "  Data saved in file '" << data_filename << "'\n";
  //
  //  Create the command file.
  //
  command_unit.open(command_filename.c_str());
  command_unit << "set term png\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set size ratio " << ratio << "\n";
  command_unit << "set timestamp\n";
  command_unit << "unset key\n";
  command_unit << "set output 'symmetric.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set title 'NACA Symmetric Airfoil'\n";
  command_unit << "plot '" << data_filename << "' using 1:2 with lines lw 3\n";
  command_unit << "quit\n";

  command_unit.close();

  std::cout << "  Created command file '" << command_filename << "'\n";

  return;
}

void test02() {
  //****************************************************************************80
  //
  //  Purpose:
  //
  //    TEST02 tests NACA4_CAMBERED.
  //
  //  Licensing:
  //
  //    This code is distributed under the MIT license.
  //
  //  Modified:
  //
  //    21 May 2014
  //
  //  Author:
  //
  //    John Burkardt
  //

  std::string command_filename = "cambered_commands.txt";
  std::ofstream command_unit;
  std::string data_filename = "cambered_data.txt";
  int n = 51;
  std::vector<double> xl;
  std::vector<double> xu;
  std::vector<double> xy;
  std::vector<double> yl;
  std::vector<double> yu;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  NACA4_CAMBERED evaluates (xu,yu) and (xl,yl) for a NACA\n";
  std::cout << "  cambered airfoil defined by a 4-digit code.\n";

  const auto m = 0.02;
  const auto p = 0.4;
  const auto t = 0.12;
  const auto c = 10.0;

  const auto xc = r8vec_linspace_new(n, 0.0, c);

  xu.resize(n);
  xl.resize(n);
  yu.resize(n);
  yl.resize(n);

  naca4_cambered(m, p, t, c, n, xc, xu, yu, xl, yl);
  //
  //  Reorganize data into a single object.
  //
  xy.resize(2 * 2 * n);

  for (auto i = 0; i < n; i++) {
    xy[0 + i * 2] = xl[i];
    xy[1 + i * 2] = yl[i];
  }
  for (auto i = 0; i < n; i++) {
    xy[0 + (n + i) * 2] = xu[n - 1 - i];
    xy[1 + (n + i) * 2] = yu[n - 1 - i];
  }
  //
  //  Determine size ratio.
  //
  const auto x_min = std::fmin(r8vec_min(n, xl), r8vec_min(n, xu));
  const auto x_max = std::fmax(r8vec_max(n, xl), r8vec_max(n, xu));
  const auto y_min = std::fmin(r8vec_min(n, yl), r8vec_min(n, yu));
  const auto y_max = std::fmax(r8vec_max(n, yl), r8vec_max(n, yu));
  const auto ratio = (y_max - y_min) / (x_max - x_min);
  //
  //  Save data to a file.
  //
  r8mat_write(data_filename, 2, 2 * n, xy);
  std::cout << "  Data saved in file '" << data_filename << "'\n";
  //
  //  Create the command file.
  //
  command_unit.open(command_filename.c_str());
  command_unit << "set term png\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set size ratio " << ratio << "\n";
  command_unit << "set timestamp\n";
  command_unit << "unset key\n";
  command_unit << "set output 'cambered.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set title 'NACA Cambered Airfoil'\n";
  command_unit << "plot '" << data_filename << "' using 1:2 with lines lw 3\n";
  command_unit << "quit\n";

  command_unit.close();

  std::cout << "  Created command file '" << command_filename << "'\n";

  return;
}

TEST(NacaTest, SymmetricAirfoil) {
  test01();
  // Add assertions to validate the output if necessary
}

TEST(NacaTest, CamberedAirfoil) {
  test02();
  // Add assertions to validate the output if necessary
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
