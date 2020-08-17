#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#include "common.h"
#include "mean_field.h"
#include "self_consistent.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <boost/multi_array.hpp>
#include <cmath>
#include <complex.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace boost;
using namespace Eigen;

int main() {
  double u = 10.0;
  double up = 5.0;
  double j = 0.0;
  int max_kxy = 100;
  auto phys_val = scc(100, u, up, j);
  multi_array<complex<double>, 2> n = get<0>(phys_val);
  multi_array<complex<double>, 2> m = get<1>(phys_val);
  double mu = get<2>(phys_val);
  string oppath = "../data/op.txt";
  ofstream opfile(oppath);
  rep(inr, 4) {
    rep(inc, 4) { opfile << n[inr][inc] << endl; }
  }
  rep(inr, 4) {
    rep(inc, 4) { opfile << m[inr][inc] << endl; }
  }
  opfile << mu << endl;
  opfile.close();
}
