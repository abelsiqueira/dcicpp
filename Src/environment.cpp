#include "environment.h"

namespace DCI {
  using namespace base_matrices;

//  Environment::Environment (const base_common & env) : base_common (env) {
//  }

//  Environment::Environment (const Environment & env) : base_common (env) {
//  }

  Environment::Environment () : base_common () {
//    useSupernodal();
    useSimplicial();
    common->quick_return_if_not_posdef = 1;
    common->dbound = 1e-6;
    common->print = 0;
  }

  Environment::~Environment () {
  }

  Bool Environment::IsPosDef () {
    if (common->status > 0)
      return dciFalse;
    else
      return dciTrue;
  }

}
