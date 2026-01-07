#include <fstream>
#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "math/include/spline.h"

namespace feasst {

double enlj(const double r) { return 4*(std::pow(r, -12) - std::pow(r, -6)); }

TEST(LinearSpline1D, LJ) {
  const double lower = 0.75, upper = 3;
  std::vector<double> xs = range(lower, upper, 8);
  std::vector<double> x2s = range(lower+1e-8, upper-1e-8, 1e3);
  std::vector<double> ys;//, y2s;
  std::ofstream file("tmp/spline.txt");
  for (const double x : xs) ys.push_back(enlj(x));
  LinearSpline1D lsp(ys);
  LinearSpline1D lsp2 = test_serialize(lsp);
  QuadraticSpline1D qsp(ys);
  QuadraticSpline1D qsp2 = test_serialize(qsp);
  // gnuplot set yrange [-1:0]
  // gnuplot p 'spline.txt' u 1:2 w l, 'spline.txt' u 1:3 w l, 'spline.txt' u 1:4 w l
  for (const double x : x2s) {
    const double z = (x - lower) / (upper - lower);
    file << x << " " << enlj(x) << " " << lsp2.interpolate(z) << " " << qsp2.interpolate(z) << std::endl;
  }
  EXPECT_NEAR(lsp2.interpolate(0.45), -0.13699525391559644, NEAR_ZERO);
  EXPECT_NEAR(qsp2.interpolate(0.45), -0.13117527694647810, NEAR_ZERO);
}

double gauss2d(const double x, const double y) { return std::exp(-(x*x + y*y)); }

TEST(LinearSpline1D, gaussian) {
  const double xlo = -2, xhi = 2, ylo = xlo, yhi = xhi;
  std::vector<double> xs = range(xlo, xhi, 8);
  std::vector<double> ys = range(ylo, yhi, 1e1);
  std::vector<double> x2s = range(xlo+1e-8, xhi-1e-8, 8e1);
  std::vector<double> y2s = range(ylo+1e-8, yhi-1e-8, 1e2);
  std::vector<std::vector<double> > fs, fxy, fxy2;
  for (const double x : xs) {
    std::vector<double> fst;
    for (const double y : ys) {
      const double z = gauss2d(x, y);
      fst.push_back(z);
      fxy.push_back({x, y, z});
    }
    fs.push_back(fst);
  }
  for (const double x : x2s) {
    for (const double y : y2s) {
      fxy2.push_back({x, y, gauss2d(x, y)});
    }
  }
  LinearSpline2D lsp(fs);
  std::ofstream file("tmp/gauss.txt");
  for (std::vector<double> xy : fxy2) {
    const double xs = (xy[0] - xlo)/(xhi - xlo);
    const double ys = (xy[1] - ylo)/(yhi - ylo);
    file << feasst_str(xy) << "," << lsp.interpolate(xs, ys) << std::endl;
  }
}

}  // namespace feasst
