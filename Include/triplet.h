#ifndef dci_triplet_h
#define dci_triplet_h

#include "vector.h"

namespace DCI {
  using namespace base_matrices;

  class Triplet : public base_triplet {
    public:
      Triplet ();
      Triplet (const Triplet &);
      Triplet (const Environment &, size_t = 0, size_t = 0, size_t = 0, int = 0); //size and nnz
      ~Triplet ();

    private:
      friend class Interface;
      double * get_doublex () const {
        return static_cast < double * > (triplet->x);
      }

      Int * get_linti () const {
        return static_cast < Int * > (triplet->i);
      }

      Int * get_lintj () const {
        return static_cast < Int * > (triplet->j);
      }

      size_t * get_pnnz () const {
        return &triplet->nnz;
      }
  };
}

#endif
