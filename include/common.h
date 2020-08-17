#ifndef COMMON_H
#define COMMON_H
#include <Eigen/Eigenvalues>
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace boost;

double fermi(const double &band, const double &temperature);
double fermi(const double &band, const int &temperature);
Matrix<complex<double>, 8, 8> h1(const double &kx, const double &ky);
Matrix<complex<double>, 8, 8> h2(const multi_array<complex<double>, 2> &n,
                                 const multi_array<complex<double>, 2> &m,
                                 const double &u, const double up,
                                 const double j);
Matrix<complex<double>, 4, 4> hintra(const int g,
                                     const multi_array<complex<double>, 2> &n,
                                     const multi_array<complex<double>, 2> &m,
                                     const double &u, const double up,
                                     const double j);
Matrix<complex<double>, 4, 4> hinter(const int g, const int gp,
                                     const multi_array<complex<double>, 2> &n,
                                     const multi_array<complex<double>, 2> &m,
                                     const double &u, const double up,
                                     const double j);
// template <typename T> class EPT {
// public:
//  double E;
//  T P;
//};
//
// template <typename T>
// T tetrahedron(const double &mu, multi_array<T, 2> physval,
//              multi_array<double, 2> Energys_temp, const int &max_kxy) {
//  vector<T> physval_vec(2 * max_kxy * max_kxy);
//  rep(triangle_ud, 2) {
//#pragma omp parallel for
//    rep(ixy, max_kxy * max_kxy) {
//      int ix = ixy / max_kxy;
//      int iy = ixy % max_kxy;
//      vector<EPT<T>> EP_k(3);
//      T physval_temp = 0.0;
//      if (triangle_ud == 0) {
//        EP_k[0].E = Energys_temp[ix][iy];
//        EP_k[1].E = Energys_temp[ix + 1][iy];
//        EP_k[2].E = Energys_temp[ix][iy + 1];
//        EP_k[0].P = physval[ix][iy];
//        EP_k[1].P = physval[ix + 1][iy];
//        EP_k[2].P = physval[ix][iy + 1];
//      } else {
//        EP_k[0].E = Energys_temp[ix + 1][iy + 1];
//        EP_k[1].E = Energys_temp[ix + 1][iy];
//        EP_k[2].E = Energys_temp[ix][iy + 1];
//        EP_k[0].P = physval[ix + 1][iy + 1];
//        EP_k[1].P = physval[ix + 1][iy];
//        EP_k[2].P = physval[ix][iy + 1];
//      }
//      sort(EP_k.begin(), EP_k.end(),
//           [](const EPT<T> &x, const EPT<T> &y) { return x.E < y.E; });
//      double S;
//      if (mu < EP_k[0].E) {
//        continue;
//      } else if (mu < EP_k[1].E) {
//        S = 0.5 * norm(mu - EP_k[0].E) / (EP_k[1].E - EP_k[0].E) /
//            (EP_k[2].E - EP_k[0].E);
//        T PM1 = EP_k[0].P + (EP_k[1].P - EP_k[0].P) * (mu - EP_k[0].E) /
//                                (EP_k[1].E - EP_k[0].E);
//        T PM2 = EP_k[0].P + (EP_k[2].P - EP_k[0].P) * (mu - EP_k[0].E) /
//                                (EP_k[2].E - EP_k[0].E);
//        physval_temp += S / 3.0 * (EP_k[0].P + PM1 + PM2);
//      } else if (mu < EP_k[2].E) {
//        T PM0 = EP_k[0].P + (EP_k[2].P - EP_k[0].P) * (mu - EP_k[0].E) /
//                                (EP_k[2].E - EP_k[0].E);
//        T PM1 = EP_k[1].P + (EP_k[2].P - EP_k[1].P) * (mu - EP_k[1].E) /
//                                (EP_k[2].E - EP_k[1].E);
//        S = 0.5 * norm(EP_k[2].E - mu) / (EP_k[2].E - EP_k[0].E) /
//            (EP_k[2].E - EP_k[1].E);
//        physval_temp += 0.5 * (EP_k[0].P + EP_k[1].P + EP_k[2].P) / 3.0 -
//                        S / 3.0 * (EP_k[2].P + PM0 + PM1);
//      } else if (EP_k[2].E < mu) {
//        physval_temp += 0.5 * (EP_k[0].P + EP_k[1].P + EP_k[2].P) / 3.0;
//      }
//      physval_vec[max_kxy * max_kxy * triangle_ud + max_kxy * ix + iy] =
//          physval_temp;
//    }
//  }
//  T phys_val_zero = 0.0;
//  return 2.0 *
//         accumulate(physval_vec.begin(), physval_vec.end(), phys_val_zero) /
//         (1.0 * max_kxy * max_kxy);
//}
template <typename T> class EPT {
public:
  double E;
  T P;
};

complex<double> tetrahedron(const double &mu,
                            multi_array<complex<double>, 2> physval,
                            multi_array<double, 2> Energys_temp,
                            const int &max_kxy);
#endif
