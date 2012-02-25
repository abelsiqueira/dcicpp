#include "environment.h"

namespace DCI {
  using namespace base_matrices;

  Bool Environment::IsPosDef () {
    if (common->status == CHOLMOD_NOT_POSDEF)
      return dciFalse;
    else
      return dciTrue;
  }

//  Environment::Environment (const base_common & env) : base_common (env) {
//  }

//  Environment::Environment (const Environment & env) : base_common (env) {
//  }

  Environment::Environment () : base_common () {
  }

  Environment::~Environment () {
  }
}
