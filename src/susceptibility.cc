#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#include "susceptibility.h"
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
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;
using namespace boost;
using namespace Eigen;

complex<double> susceptibility::calc_susceptibility(
    const double &temprature, const double &mu, double omega, const int &sigma1,
    const int &ikappa, const int &sigma2, const int &ilambda, const int &l1,
    const int &sigma3, const int &imu, const int &sigma4, const int &inu,
    const int &l2) {
  vector<complex<double>> result(4 * 64 * (max_kxy) * (max_kxy));
  complex<double> j = 1.0i;
#pragma omp parallel for
  rep(i, (max_kxy) * (max_kxy)) {
    int ix = i / (max_kxy);
    int iy = i % (max_kxy);
    rep(mq, 2) {
      rep(nq, 2) {
        int m1 = (mq) % 2;
        int m2 = (mq + l1) % 2;
        int m3 = (mq + nq + l2) % 2;
        int m4 = (mq + nq) % 2;
        int i1 = sigma1 + 2 * m1 + 4 * ikappa;
        int i2 = sigma2 + 2 * m2 + 4 * ilambda;
        int i3 = sigma3 + 2 * m3 + 4 * imu;
        int i4 = sigma4 + 2 * m4 + 4 * inu;
        rep(epsilon1, 8) {
          rep(epsilon2, 8) {
            result[4 * 64 * i + 2 * 64 * mq + 64 * nq + 8 * epsilon1 +
                   epsilon2] =
                -(fermi(energyk[0][epsilon1][ix][iy] - mu, 0) -
                  fermi(energyk[1][epsilon2][ix][iy] - mu, 0)) *
                eigenvectork[0][i2][epsilon1][ix][iy] *
                conj(eigenvectork[0][i3][epsilon1][ix][iy]) *
                eigenvectork[1][i4][epsilon2][ix][iy] *
                conj(eigenvectork[1][i1][epsilon2][ix][iy]) /
                (energyk[0][epsilon1][ix][iy] - energyk[1][epsilon2][ix][iy] -
                 (omega + 0.01 * j));
          }
        }
      }
    }
  }
  complex<double> czero = 0.0;
  return accumulate(result.begin(), result.end(), czero) /
         (2.0 * max_kxy * max_kxy);
}

void susceptibility::calc_energy(double u, double up, double j,
                                 multi_array<complex<double>, 2> n,
                                 multi_array<complex<double>, 2> m) {
  Matrix<complex<double>, 8, 8> hi = h2(n, m, u, up, j);
#pragma omp parallel for
  rep(i, (max_kxy) * (max_kxy)) {
    int i0 = i / (max_kxy);
    int i1 = i % (max_kxy);
    double kx = M_PI * (i0 - i1) / max_kxy;
    double ky = -M_PI + M_PI * (i0 + i1) / max_kxy;
    Matrix<complex<double>, 8, 8> hamiltonian0 = h1(kx + q, ky + q) + hi;
    Matrix<complex<double>, 8, 8> hamiltonian1 = h1(kx, ky) + hi;
    SelfAdjointEigenSolver<Matrix<complex<double>, 8, 8>> es0(hamiltonian0);
    SelfAdjointEigenSolver<Matrix<complex<double>, 8, 8>> es1(hamiltonian1);
    rep(ei, 8) {
      energyk[0][ei][i0][i1] = es0.eigenvalues()(ei);
      energyk[1][ei][i0][i1] = es1.eigenvalues()(ei);
      rep(gi, 8) {
        eigenvectork[0][gi][ei][i0][i1] = es0.eigenvectors()(gi, ei);
        eigenvectork[1][gi][ei][i0][i1] = es1.eigenvectors()(gi, ei);
        /* eigenvectork[0][gi][ei][i0][i1] = es0.eigenvectors()(ei, gi); */
        /* eigenvectork[1][gi][ei][i0][i1] = es1.eigenvectors()(ei, gi); */
      }
    }
  }
}

/*
index = sigma1 + 2*icf1 + 4*sigma2 + 8* icf2 + 16*l1
*/

vector<int> i2l32(int index) {
  int ind = index;
  int l1 = ind / 16;
  ind = ind % 16;
  int icf2 = ind / 8;
  ind = ind % 8;
  int sigma2 = ind / 4;
  ind = ind % 4;
  int icf1 = ind / 2;
  ind = ind % 2;
  int sigma1 = ind;
  vector<int> result = {sigma1, icf1, sigma2, icf2, l1};
  return result;
}

/*
index = icf1 + 2*icf2 + 4*l1
*/
vector<int> i2l8(int index) {
  int ind = index;
  int l1 = ind / 4;
  ind = ind % 4;
  int icf2 = ind / 2;
  ind = ind % 2;
  int icf1 = ind;
  vector<int> result = {icf1, icf2, l1};
  return result;
}

int labels2index(int sigma1, int icf1, int sigma2, int icf2, int l1) {
  return sigma1 + 2 * icf1 + 4 * sigma2 + 8 * icf2 + 16 * l1;
}

double v(int a, int b, int c, int d, double u, double up, double j) {
  double result;
  if (a == b and b == c and c == d) {
    result = u;
  } else if (a == c and c != b and b == d) {
    result = j;
  } else if (a == b and b != c and c == d) {
    result = j;
  } else if (a == d and d != b and b == c) {
    result = up;
  } else {
    result = 0.0;
  }
  return result;
}

int main() {
  double u = 10.0;
  double up = 5.0;
  double j = 0.0;
  int max_kxy = 100;
  multi_array<complex<double>, 2> n, m;
  n.resize(extents[4][4]);
  m.resize(extents[4][4]);
  double mu;
  string oppath = "../data/op.txt";
  ifstream opfile(oppath);
  rep(inr, 4) {
    rep(inc, 4) { opfile >> n[inr][inc]; }
  }
  rep(inr, 4) {
    rep(inc, 4) { opfile >> m[inr][inc]; }
  }
  opfile >> mu;
  opfile.close();
  string dirname = "brydon3";
  string chi0fname = "../data/" + dirname + "/chi0.txt";
  string chifname = "../data/" + dirname + "/chi.txt";
  ofstream chi0file(chi0fname);
  ofstream chifile(chifname);
  chi0file << "q,omega,chi0" << endl;
  chifile << "q,omega,chi" << endl;
  // multi_array<complex<double>, 2> chi(extents[2][2]);
  /* Matrix<complex<double>, 8, 8> chi0; */
  MatrixXcd chi0 = MatrixXcd::Zero(8, 8);
  const int max_iomega = 40;
  const int max_iq = 20;
  rep(iq, max_iq + 1) {
    double q = M_PI * iq / max_iq;
    susceptibility susceptibility_solver(q, max_kxy);
    susceptibility_solver.calc_energy(u, up, j, n, m);
    rep(iomega, max_iomega + 1) {
      double omega = 4.0 * iomega / max_iomega;
      rep(irow, 8) {
        vector<int> rl = i2l8(irow);
        rep(icol, 8) {
          vector<int> cl = i2l8(icol);
          chi0(irow, icol) = susceptibility_solver.calc_susceptibility(
              0, mu, omega, 1, rl[0], 0, rl[1], rl[2], 0, cl[0], 1, cl[1],
              cl[2]);
        }
      }
      cout << "##" << endl;
      cout << "q = " << q << " omega = " << omega << endl;
      rep(irow, 8) {
        rep(icol, 8) {
          chi0file << q << "," << omega << "," << real(chi0(irow, icol)) << endl
                   << q << "," << omega << "," << imag(chi0(irow, icol))
                   << endl;
        }
      }
      /*
      index = sigma1 + 2*icf1 + 4*sigma2 + 8* icf2 + 16*l1
      */
      /* Matrix<complex<double>, 8, 8> dysonmat; */
      MatrixXcd dysonmat = MatrixXcd::Zero(8, 8);
      rep(irow, 8) {
        dysonmat(irow, irow) = 1.0;
        rep(icol, 8) {
          vector<int> lc = i2l8(icol);
          rep(ialpha, 2) {
            rep(ibeta, 2) {
              dysonmat(irow, icol) +=
                  -chi0(irow, ialpha + 2 * ibeta + 4 * lc[2]) *
                  v(ialpha, ibeta, lc[0], lc[1], u, up, j);
            }
          }
        }
      }
      Matrix<complex<double>, 8, 8> chi = dysonmat.inverse() * chi0;
      rep(irow, 8) {
        rep(icol, 8) {
          chifile << q << "," << omega << "," << real(chi(irow, icol)) << endl
                  << q << "," << omega << "," << imag(chi(irow, icol)) << endl;
        }
      }
    }
  }
  chi0file.close();
  chifile.close();
}
