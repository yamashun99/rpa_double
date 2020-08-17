#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#include "common.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

double fermi(const double &band, const double &temperature) {
  if (temperature > 0) {
    return 1 / (exp(band / temperature) + 1);
  } else if (temperature == 0) {
    if (band < 0) {
      return 1;
    } else if (band > 0) {
      return 0;
    } else {
      return 0.5;
    }
  }
  return -1;
}

double fermi(const double &band, const int &temperature) {
  if (band < 0) {
    return 1;
  } else if (band > 0) {
    return 0;
  } else {
    return 0.5;
  }
  return -1;
}

Matrix<complex<double>, 8, 8> h1(const double &kx, const double &ky) {
  double t = 1.0;
  double d = 6.0;
  Matrix<complex<double>, 8, 8> result;
  result(0, 0) = 2.0 * t * (cos(kx) + cos(ky)) - d / 2.0;
  result(2, 2) = -2.0 * t * (cos(kx) + cos(ky)) - d / 2.0;
  result(4, 4) = 2.0 * t * (cos(kx) + cos(ky)) + d / 2.0;
  result(6, 6) = -2.0 * t * (cos(kx) + cos(ky)) + d / 2.0;
  rep(ig, 4) { result(2 * ig + 1, 2 * ig + 1) = result(2 * ig, 2 * ig); }
  return result;
}

Matrix<complex<double>, 8, 8> h2(const multi_array<complex<double>, 2> &n,
                                 const multi_array<complex<double>, 2> &m,
                                 const double &u, const double up,
                                 const double j) {
  Matrix<complex<double>, 8, 8> result;
  rep(g, 2) {
    rep(ri, 4) {
      rep(cj, 4) {
        result(4 * g + ri, 4 * g + cj) = hintra(g, n, m, u, up, j)(ri, cj);
      }
    }
  }
  rep(g, 2) {
    rep(gp, 2) {
      if (g <= gp) {
        continue;
      }
      rep(ri, 4) {
        rep(cj, 4) {
          result(4 * g + ri, 4 * gp + cj) =
              hinter(g, gp, n, m, u, up, j)(ri, cj);
        }
      }
    }
  }
  return result;
}

Matrix<complex<double>, 4, 4> hintra(const int g,
                                     const multi_array<complex<double>, 2> &n,
                                     const multi_array<complex<double>, 2> &m,
                                     const double &u, const double up,
                                     const double j) {
  vector<complex<double>> rc(6);
  rc[0] = u * n[2 * g + 1][2 * g + 1];
  rep(gp, 2) {
    if (gp != g) {
      rc[0] += up * n[2 * gp + 1][2 * gp + 1] + (up - j) * n[2 * gp][2 * gp];
    }
  }
  rc[1] = -u * n[2 * g][2 * g + 1];
  rep(gp, 2) {
    if (gp != g) {
      rc[1] += -j * n[2 * gp][2 * gp + 1];
    }
  }
  rc[2] = u * n[2 * g][2 * g];
  rep(gp, 2) {
    if (gp != g) {
      rc[2] += up * n[2 * gp][2 * gp] + (up - j) * n[2 * gp + 1][2 * gp + 1];
    }
  }
  rc[3] = u * m[2 * g + 1][2 * g + 1];
  rep(gp, 2) {
    if (gp != g) {
      rc[3] += up * m[2 * gp + 1][2 * gp + 1] + (up - j) * m[2 * gp][2 * gp];
    }
  }
  rc[4] = -u * m[2 * g][2 * g + 1];
  rep(gp, 2) {
    if (gp != g) {
      rc[4] += -j * m[2 * gp][2 * gp + 1];
    }
  }
  rc[5] = u * m[2 * g][2 * g];
  rep(gp, 2) {
    if (gp != g) {
      rc[5] += up * m[2 * gp][2 * gp] + (up - j) * m[2 * gp + 1][2 * gp + 1];
    }
  }
  Matrix<complex<double>, 4, 4> result;
  // clang-format off
   result << rc[0], 0,           0,     0,
             rc[1], rc[2],       0,     0,
             rc[3], conj(rc[4]), rc[0], 0,
             rc[4], rc[5],       rc[1], rc[2];
  // clang-format on
  return result;
}

Matrix<complex<double>, 4, 4> hinter(const int g, const int gp,
                                     const multi_array<complex<double>, 2> &n,
                                     const multi_array<complex<double>, 2> &m,
                                     const double &u, const double up,
                                     const double j) {
  vector<complex<double>> rc(8);
  rc[0] = j * n[2 * gp + 1][2 * g + 1] + j * n[2 * g + 1][2 * gp + 1] -
          (up - j) * n[2 * gp][2 * g];
  rc[1] = -j * n[2 * g + 1][2 * gp] - up * n[2 * gp + 1][2 * g];
  rc[2] = -j * n[2 * g][2 * gp + 1] - up * n[2 * gp][2 * g + 1];
  rc[3] = j * n[2 * gp][2 * g] + j * n[2 * g][2 * gp] -
          (up - j) * n[2 * gp + 1][2 * g + 1];
  rc[4] = j * m[2 * gp + 1][2 * g + 1] + j * m[2 * g + 1][2 * gp + 1] -
          (up - j) * m[2 * gp][2 * g];
  rc[5] = -j * m[2 * g + 1][2 * gp] - up * m[2 * gp + 1][2 * g];
  rc[6] = -j * m[2 * g][2 * gp + 1] - up * m[2 * gp][2 * g + 1];
  rc[7] = j * m[2 * gp][2 * g] + j * m[2 * g][2 * gp] -
          (up - j) * m[2 * gp + 1][2 * g + 1];
  Matrix<complex<double>, 4, 4> result;
  // clang-format off
   result << rc[0], rc[1], rc[4], rc[5],
             rc[2], rc[3], rc[6], rc[7],
             rc[4], rc[5], rc[0], rc[1],
             rc[6], rc[7], rc[2], rc[3];
  // clang-format on
  return result;
}

complex<double> tetrahedron(const double &mu,
                            multi_array<complex<double>, 2> physval,
                            multi_array<double, 2> Energys_temp,
                            const int &max_kxy) {
  vector<complex<double>> physval_vec(2 * max_kxy * max_kxy);
  rep(triangle_ud, 2) {
#pragma omp parallel for
    rep(ixy, max_kxy * max_kxy) {
      int ix = ixy / max_kxy;
      int iy = ixy % max_kxy;
      vector<EPT<complex<double>>> EP_k(3);
      complex<double> physval_temp = 0.0;
      if (triangle_ud == 0) {
        EP_k[0].E = Energys_temp[ix][iy];
        EP_k[1].E = Energys_temp[ix + 1][iy];
        EP_k[2].E = Energys_temp[ix][iy + 1];
        EP_k[0].P = physval[ix][iy];
        EP_k[1].P = physval[ix + 1][iy];
        EP_k[2].P = physval[ix][iy + 1];
      } else {
        EP_k[0].E = Energys_temp[ix + 1][iy + 1];
        EP_k[1].E = Energys_temp[ix + 1][iy];
        EP_k[2].E = Energys_temp[ix][iy + 1];
        EP_k[0].P = physval[ix + 1][iy + 1];
        EP_k[1].P = physval[ix + 1][iy];
        EP_k[2].P = physval[ix][iy + 1];
      }
      sort(EP_k.begin(), EP_k.end(),
           [](const EPT<complex<double>> &x, const EPT<complex<double>> &y) {
             return x.E < y.E;
           });
      double S;
      if (mu < EP_k[0].E) {
        continue;
      } else if (mu < EP_k[1].E) {
        S = 0.5 * norm(mu - EP_k[0].E) / (EP_k[1].E - EP_k[0].E) /
            (EP_k[2].E - EP_k[0].E);
        complex<double> PM1 = EP_k[0].P + (EP_k[1].P - EP_k[0].P) *
                                              (mu - EP_k[0].E) /
                                              (EP_k[1].E - EP_k[0].E);
        complex<double> PM2 = EP_k[0].P + (EP_k[2].P - EP_k[0].P) *
                                              (mu - EP_k[0].E) /
                                              (EP_k[2].E - EP_k[0].E);
        physval_temp += S / 3.0 * (EP_k[0].P + PM1 + PM2);
      } else if (mu < EP_k[2].E) {
        complex<double> PM0 = EP_k[0].P + (EP_k[2].P - EP_k[0].P) *
                                              (mu - EP_k[0].E) /
                                              (EP_k[2].E - EP_k[0].E);
        complex<double> PM1 = EP_k[1].P + (EP_k[2].P - EP_k[1].P) *
                                              (mu - EP_k[1].E) /
                                              (EP_k[2].E - EP_k[1].E);
        S = 0.5 * norm(EP_k[2].E - mu) / (EP_k[2].E - EP_k[0].E) /
            (EP_k[2].E - EP_k[1].E);
        physval_temp += 0.5 * (EP_k[0].P + EP_k[1].P + EP_k[2].P) / 3.0 -
                        S / 3.0 * (EP_k[2].P + PM0 + PM1);
      } else if (EP_k[2].E < mu) {
        physval_temp += 0.5 * (EP_k[0].P + EP_k[1].P + EP_k[2].P) / 3.0;
      }
      physval_vec[max_kxy * max_kxy * triangle_ud + max_kxy * ix + iy] =
          physval_temp;
    }
  }
  complex<double> phys_val_zero = 0.0;
  return accumulate(physval_vec.begin(), physval_vec.end(), phys_val_zero) /
         (2.0 * max_kxy * max_kxy);
}
