#ifndef MEAN_FIELD_H
#define MEAN_FIELD_H
#include <Eigen/Eigenvalues>
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>
using namespace std;
using namespace Eigen;
using namespace boost;

class mf {
private:
  int max_kxy;
  boost::multi_array<complex<double>, 5> nk;
  boost::multi_array<complex<double>, 5> mk;

public:
  boost::multi_array<double, 3> energyk;
  multi_array<complex<double>, 4> eigenvectork;
  multi_array<complex<double>, 3>
  main_mf_loop(const multi_array<complex<double>, 2> &n_guess,
               const multi_array<complex<double>, 2> &m_guess, const double &u,
               const double &up, const double &j, const double &mu);
  void calc_phys_k(multi_array<complex<double>, 2> n,
                   multi_array<complex<double>, 2> m, double u, double up,
                   double j);
  double pn(double mu);
  double cp(double n, int ind);
  mf(const int &max_kxy) {
    this->max_kxy = max_kxy;
    nk.resize(extents[4][4][8][max_kxy + 1][max_kxy + 1]);
    mk.resize(extents[4][4][8][max_kxy + 1][max_kxy + 1]);
    energyk.resize(extents[8][max_kxy + 1][max_kxy + 1]);
    eigenvectork.resize(extents[8][8][max_kxy + 1][max_kxy + 1]);
  }
};
#endif
