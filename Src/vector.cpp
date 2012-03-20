#include "vector.h"

namespace DCI {
  using namespace base_matrices;

  Vector::Vector (const base_dense & v) : base_dense (v) {
  }

  Vector::Vector (const Vector & v) : base_dense (v) {
  }

  Vector::Vector (const Environment & env, size_t size, double * v) : base_dense (env, size, 1) {
    if ( (v != 0) && (size > 0) ) {
      double * px = static_cast < double * > (dense->x);

      for (size_t i = 0; i < size; i++)
        px[i] = v[i];
    }
  }

  Vector::Vector (const Environment & env, const std::vector < double > & V) : base_dense (env, V.size(), 1) {
    size_t n = V.size();
    double * px = static_cast < double * > (dense->x);

    for (size_t i = 0; i < n; i++)
      px[i] = V[i];

  }

  Vector::~Vector () {
  }

  void Vector::reset (size_t n, double a) {
    // If it is allocated but the size is not n, free
    if ( (dense != 0) && (dense->nrow != n) ) {
      CHOLMOD(free_dense) (&dense, get_cholmod_common());
    }
    // If freed before, or never allocated
    if ( dense == 0 )
      dense = CHOLMOD(allocate_dense) (n, 1, n, CHOLMOD_REAL, get_cholmod_common());
    // Now, dense has size n
    double * px = static_cast < double * > (dense->x);
    for (size_t i = 0; i < n; i++)
      px[i] = a;

  }

  //Access
  double Vector::get (size_t pos) {
    if ( (pos == 0) || (pos > dense->nrow) )
      return 0;
    double * px = static_cast < double * > (dense->x);
    return px[pos - 1];
  }

  void Vector::set (size_t pos, double value) {
    if ( (pos == 0) || (pos > dense->nrow) )
      return;
    double * px = static_cast < double * > (dense->x);
    px[pos - 1] = value;
  }

  double Vector::operator() (size_t pos) {
    if ( (pos == 0) || (pos > dense->nrow) )
      return 0;
    double * px = static_cast < double * > (dense->x);
    return px[pos - 1];
  }

}
