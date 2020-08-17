#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define EIGEN_NO_DEBUG
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_MPL2_ONLY
#include "mean_field.h"
#include "common.h"
#include <boost/multi_array.hpp>
#include <cxx-prettyprint/prettyprint.hpp>
#include <iostream>
#include <numeric>

using namespace std;
using namespace Eigen;
using namespace boost;

multi_array<complex<double>, 3>
mf::main_mf_loop(const multi_array<complex<double>, 2> &n_guess,
                 const multi_array<complex<double>, 2> &m_guess,
                 const double &u, const double &up, const double &j,
                 const double &mu) {
  multi_array<complex<double>, 3> result;
  result.resize(extents[2][4][4]);
  calc_phys_k(n_guess, m_guess, u, up, j);
  multi_array<complex<double>, 2> n;
  multi_array<complex<double>, 2> m;
  n.resize(extents[4][4]);
  m.resize(extents[4][4]);
  rep(ei, 8) {
    rep(g1, 4) {
      rep(g2, 4) {
        n[g1][g2] += tetrahedron(mu, nk[g1][g2][ei], energyk[ei], max_kxy);
        m[g1][g2] += tetrahedron(mu, mk[g1][g2][ei], energyk[ei], max_kxy);
      }
    }
  }
  result[0] = n;
  result[1] = m;
  return result;
}

void mf::calc_phys_k(multi_array<complex<double>, 2> n,
                     multi_array<complex<double>, 2> m, double u, double up,
                     double j) {
  Matrix<complex<double>, 8, 8> hi = h2(n, m, u, up, j);
#pragma omp parallel for
  rep(i, (max_kxy + 1) * (max_kxy + 1)) {
    int i0 = i / (max_kxy + 1);
    int i1 = i % (max_kxy + 1);
    double kx = M_PI * (i0 - i1) / max_kxy;
    double ky = -M_PI + M_PI * (i0 + i1) / max_kxy;
    Matrix<complex<double>, 8, 8> hamiltonian = h1(kx, ky) + hi;
    SelfAdjointEigenSolver<Matrix<complex<double>, 8, 8>> es(hamiltonian);
    rep(ei, 8) {
      energyk[ei][i0][i1] = es.eigenvalues()(ei);
      rep(g1, 2) {
        rep(g2, 2) {
          rep(s1, 2) {
            rep(s2, 2) {
              nk[2 * g1 + s1][2 * g2 + s2][ei][i0][i1] =
                  1.0 * (conj(es.eigenvectors()(4 * g1 + s1, ei)) *
                             es.eigenvectors()(4 * g2 + s2, ei) +
                         conj(es.eigenvectors()(4 * g1 + s1 + 2, ei)) *
                             es.eigenvectors()(4 * g2 + s2 + 2, ei));
              mk[2 * g1 + s1][2 * g2 + s2][ei][i0][i1] =
                  1.0 * (conj(es.eigenvectors()(4 * g1 + s1 + 2, ei)) *
                             es.eigenvectors()(4 * g2 + s2, ei) +
                         conj(es.eigenvectors()(4 * g1 + s1, ei)) *
                             es.eigenvectors()(4 * g2 + s2 + 2, ei));
            }
          }
        }
      }
    }
  }
}

double mf::cp(double inN, int ind) {
  vector<double> mu_search;
  mu_search = {-40.0, 40.0};
  vector<double> N_search = {0.0, 16.0};
  int inear;
  if (abs(inN - N_search[0]) < abs(inN - N_search[1])) {
    inear = 0;
  } else {
    inear = 1;
  }
  double epsilon = 1e-8;
  while (abs(inN - N_search[inear]) > epsilon) {
    // cout << mu_search << endl;
    // cout << N_search[0] - inN << "," << N_search[1] - inN << endl;
    sort(mu_search.begin(), mu_search.end());
    sort(N_search.begin(), N_search.end());
    //    double mu_New = mu_search[0] + (mu_search[1] - mu_search[0]) *
    //                                       (inN - N_search[0]) /
    //                                       (N_search[1] - N_search[0]);
    double mu_New = 0.5 * mu_search[0] + 0.5 * mu_search[1];
    double N_New = pn(mu_New);
    if (N_New > inN) {
      mu_search[1] = mu_New;
      N_search[1] = N_New;
    } else {
      mu_search[0] = mu_New;
      N_search[0] = N_New;
    }
    if (abs(inN - N_search[0]) < abs(inN - N_search[1])) {
      inear = 0;
    } else {
      inear = 1;
    }
  }
  return mu_search[inear];
}

double mf::pn(double mu) {
  int nb = 8;
  vector<double> N_vec(2 * nb * max_kxy * max_kxy);
  rep(triangle_ud, 2) {
    rep(l, nb) {
#pragma omp parallel for
      rep(ixy, max_kxy * max_kxy) {
        int ix = ixy / max_kxy;
        int iy = ixy % max_kxy;
        vector<double> ET(3);
        if (triangle_ud == 0) {
          ET[0] = energyk[l][ix][iy];
          ET[1] = energyk[l][ix + 1][iy];
          ET[2] = energyk[l][ix][iy + 1];
        } else {
          ET[0] = energyk[l][ix + 1][iy + 1];
          ET[1] = energyk[l][ix + 1][iy];
          ET[2] = energyk[l][ix][iy + 1];
        };
        sort(ET.begin(), ET.end());
        double N_vec_tel = 0.0;
        if (mu < ET[0]) {
          N_vec_tel = 0.0;
        } else if (mu <= ET[1]) {
          N_vec_tel =
              0.5 * norm(mu - ET[0]) / ((ET[1] - ET[0]) * (ET[2] - ET[0]));
        } else if (mu <= ET[2]) {
          N_vec_tel = 0.5 - 0.5 * norm(ET[2] - mu) /
                                ((ET[2] - ET[0]) * (ET[2] - ET[1]));
        } else if (ET[2] < mu) {
          N_vec_tel = 0.5;
        }
        N_vec[2 * max_kxy * max_kxy * l + triangle_ud * max_kxy * max_kxy +
              ixy] = N_vec_tel;
      }
    }
  }
  return 2.0 * accumulate(N_vec.begin(), N_vec.end(), 0.0) /
         (1.0 * max_kxy * max_kxy);
}
