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
  double u = 6.0;
  int max_iup = 50;
  ofstream f("../data/up_op.csv");
  f << "up,op" << endl;
  rep(iup, max_iup) {
    double up = 0.2 + 5.0 / max_iup * iup;
    double j = 0.0;
    auto phys_val = scc(100, u, up, j);
    multi_array<complex<double>, 2> n = get<0>(phys_val);
    multi_array<complex<double>, 2> m = get<1>(phys_val);
    double mu = get<2>(phys_val);
    f << up << "," << real(m[0][2]) << endl;
  }
  f.close();
}
