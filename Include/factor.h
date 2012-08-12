#ifndef dci_factor_h
#define dci_factor_h

#include "sparse.h"

namespace DCI {
  using namespace base_matrices;

//  typedef base_factor Factor;
  class Factor : public base_factor {
    public:
      Factor ();
      Factor (const Factor &);
      Factor (const Environment &, size_t = 0);
      ~Factor ();
    private:
      friend class Interface;
      double * get_doublex () const {
        return static_cast < double * > (factor->x);
      }
      Int * get_rowi () const {
        return static_cast < Int * > (factor->i);
      }
      Int * get_colp () const {
        return static_cast < Int * > (factor->p);
      }
      Int * get_colnz () const {
        return static_cast < Int * > (factor->nz);
      }
      size_t get_nnz () const {
        return factor->nzmax;
      }
      size_t get_ncol () const {
        return factor->n;
      }
  };
}

#endif
