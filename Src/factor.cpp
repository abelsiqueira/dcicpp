#include "factor.h"

namespace DCI {
  using namespace base_matrices;

  Factor::Factor (const Factor & v) : base_factor (v) {
  }

  Factor::Factor (const Environment & env, size_t n) : base_factor (env, n) {
  }

  Factor::~Factor () {
  }

}

