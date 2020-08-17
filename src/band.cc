#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#include "band.h"
#include "common.h"
#include "mean_field.h"
#include "self_consistent.h"
#include <boost/multi_array.hpp>
#include <cmath>
#include <complex>
#include <cxx-prettyprint/prettyprint.hpp>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost;

int main() {
  double u = 10.0;
  double up = 5.0;
  double j = 0.0;
  auto phys_val = scc(100, u, up, j);
  multi_array<complex<double>, 2> n = get<0>(phys_val);
  multi_array<complex<double>, 2> m = get<1>(phys_val);
  double mu = get<2>(phys_val);
  vector<vector<double>> ks(3, vector<double>(3 * 100));
  rep(i, 100) {
    double s2 = sqrt(2);
    double k = M_PI * (1.0 * i / 100);
    ks[0][i] = k;
    ks[1][i] = 0.0;
    ks[2][i] = k;
    ks[0][100 + i] = M_PI;
    ks[1][100 + i] = k;
    ks[2][100 + i] = M_PI + k;
    ks[0][2 * 100 + i] = M_PI - k;
    ks[1][2 * 100 + i] = M_PI - k;
    ks[2][2 * 100 + i] = 2.0 * M_PI + sqrt(2) * k;
  }
  ofstream f("../data/data.csv");
  int nb = 8;
  f << "k,";
  rep(i, nb) { f << "E" + to_string(i) + ","; }
  f << endl;
  rep(i, 3 * 100) {
    SelfAdjointEigenSolver<Matrix<complex<double>, Dynamic, Dynamic>> s(
        h1(ks[0][i], ks[1][i]) + h2(n, m, u, up, j));
    f << ks[2][i] << ",";
    rep(j, nb) { f << s.eigenvalues()(j) - mu << ","; }
    f << endl;
  }
  f.close();
}
