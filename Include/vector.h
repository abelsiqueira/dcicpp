#ifndef dci_vector_h
#define dci_vector_h

#include "factor.h"
#include <vector>

namespace DCI {
  using namespace base_matrices;

  class Vector : public base_dense {
    public:
      Vector ();
      Vector (const base_dense &);
      Vector (const Vector &);
      Vector (const Environment &, size_t = 0, double * = 0);
      Vector (const Environment &, const std::vector < double > &);
      ~Vector ();

      void operator= (const Vector &);

      //Access functions
      double get (size_t);
      void set (size_t, double);
      double operator() (size_t);

      //Information
      size_t size () const { return dense->nrow; };

      friend class Interface;
    private:
      double * get_doublex () const {
        return static_cast < double * > (dense->x);
      }
  };
}

#endif
