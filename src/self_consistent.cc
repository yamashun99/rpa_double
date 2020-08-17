#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#include "self_consistent.h"
#include "common.h"
#include "mean_field.h"
#include <iostream>
#include <numeric>

using namespace std;
using namespace boost;

tuple<multi_array<complex<double>, 2>, multi_array<complex<double>, 2>, double>
scc(const int &max_kxy, const double &u, const double &up, const double &j) {
  int iter_max = 10000;
  double epsilon = 1e-8;
  multi_array<complex<double>, 2> n_guess;
  multi_array<complex<double>, 2> m_guess;
  n_guess.resize(extents[4][4]);
  m_guess.resize(extents[4][4]);
  rep(ir, 4) {
    rep(ic, 4) {
      n_guess[ir][ic] = 1.0;
      m_guess[ir][ic] = 1.0;
    }
  }
  double mu;
  mf mf_loop(max_kxy);
  rep(ind, iter_max) {
    mf_loop.calc_phys_k(n_guess, m_guess, u, up, j);
    mu = mf_loop.cp(8.0, ind);
    multi_array<complex<double>, 3> nm =
        mf_loop.main_mf_loop(n_guess, m_guess, u, up, j, mu);
    multi_array<complex<double>, 2> n_guess_new = nm[0];
    multi_array<complex<double>, 2> m_guess_new = nm[1];
    bool conv = true;
    rep(ri, 4) {
      rep(ci, 4) {
        conv = conv and
               (abs(n_guess[ri][ci] - n_guess_new[ri][ci]) < epsilon) and
               (abs(n_guess[ri][ci] - n_guess_new[ri][ci]) < epsilon);
      }
    }
    if (conv) {
      break;
    }
    rep(ri, 4) {
      rep(ci, 4) {
        n_guess[ri][ci] = 0.1 * n_guess_new[ri][ci] + 0.9 * n_guess[ri][ci];
        m_guess[ri][ci] = 0.1 * m_guess_new[ri][ci] + 0.9 * m_guess[ri][ci];
      }
    }
    m_guess[0][0] = 0.0;
    m_guess[1][1] = 0.0;
    m_guess[2][2] = 0.0;
    m_guess[3][3] = 0.0;
    m_guess[0][1] = 0.0;
    m_guess[1][0] = 0.0;
    m_guess[2][3] = 0.0;
    m_guess[3][2] = 0.0;
    m_guess[0][3] = 0.0;
    m_guess[3][0] = 0.0;
    m_guess[1][2] = 0.0;
    m_guess[2][1] = 0.0;
    n_guess[0][1] = 0.0;
    n_guess[1][0] = 0.0;
    n_guess[2][3] = 0.0;
    n_guess[3][2] = 0.0;
    n_guess[0][2] = 0.0;
    n_guess[0][3] = 0.0;
    n_guess[1][2] = 0.0;
    n_guess[1][3] = 0.0;
    n_guess[2][0] = 0.0;
    n_guess[2][1] = 0.0;
    n_guess[3][0] = 0.0;
    n_guess[3][1] = 0.0;
    n_guess[1][1] = n_guess[0][0];
    n_guess[3][3] = n_guess[2][2];
    m_guess[1][3] = -m_guess[0][2];
    m_guess[3][1] = -m_guess[2][0];
    cout << "u = " << u << ", up = " << up << ", mu = " << mu << endl;
    cout << "n = " << endl;
    rep(ir, 4) {
      rep(ic, 4) { cout << n_guess[ir][ic] << " "; }
      cout << endl;
    }
    cout << "m = " << endl;
    rep(ir, 4) {
      rep(ic, 4) { cout << m_guess[ir][ic] << " "; }
      cout << endl;
    }
  }
  auto result = tie(n_guess, m_guess, mu);
  return result;
}
