#ifndef dci_environment_h
#define dci_environment_h

#include <definitions.h>
#include "base_matrices.h"

namespace DCI {
  using namespace base_matrices;

  class Vector;

  class Environment : public base_common {
    public:
      Environment ();
      Environment (const base_common &);
      Environment (const Environment &);
      ~Environment ();

      Bool IsPosDef ();
    private:
      friend class vector;
  };
}

#endif
