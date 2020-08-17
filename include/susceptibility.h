#ifndef SUSCEPTIBILITY_H
#define SUSCEPTIBILITY_H
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

using namespace std;
using namespace boost;

class susceptibility {
  double u, q;
  int max_kxy;

public:
  complex<double> calc_susceptibility(const double &temprature,
                                      const double &mu, double omega,
                                      const int &sigma1, const int &icf1,
                                      const int &sigma2, const int &icf2,
                                      const int &l1, const int &sigma3,
                                      const int &icf3, const int &sigma4,
                                      const int &icf4, const int &l2);
  void calc_energy(double u, double up, double j,
                   multi_array<complex<double>, 2> n,
                   multi_array<complex<double>, 2> m);
  multi_array<complex<double>, 5> eigenvectork;
  multi_array<double, 4> energyk;

  susceptibility(const double &q, const int &max_kxy) {
    this->max_kxy = max_kxy;
    this->q = q;
    eigenvectork.resize(extents[2][8][8][max_kxy + 1][max_kxy + 1]);
    energyk.resize(extents[2][8][max_kxy + 1][max_kxy + 1]);
  }
};

vector<int> i2l32(int index); //{sigma1, icf1, sigma2, icf2, l1}
vector<int> i2l8(int index);  //{icf1, icf2, l1}
int labels2index(int sigma1, int icf1, int sigma2, int icf2, int l1);
double v(int a, int b, int c, int d, double u, double up, double j);
#endif
