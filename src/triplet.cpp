#include "triplet.h"

namespace DCI {
  using namespace base_matrices;

  Triplet::Triplet (const Triplet & v) : base_triplet (v) {
  }

  Triplet::Triplet (const Environment & env, size_t nr, size_t nc, size_t nnz, int stype) : base_triplet (env, nr, nc, nnz, stype) {
  }

  Triplet::~Triplet () {
  }

}
