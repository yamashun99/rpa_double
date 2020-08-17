#ifndef SELF_CONSISTENT_H
#define SELF_CONSISTENT_H
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

using namespace std;
using namespace boost;

tuple<multi_array<complex<double>, 2>, multi_array<complex<double>, 2>, double>
scc(const int &max_kxy, const double &u, const double &up, const double &j);
#endif
