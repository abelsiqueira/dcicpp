#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {

  Int Interface::normalSafeguard () {
    Vector dn (*env);
    Int iout = 0;
    Bool scaleJ = scale_normal;

    call_ccfsg_xc (dciTrue, scaleJ);
    if (!is_linear)
      this->cholesky ();

    naStep (*c, dn);

#ifndef NDEBUG
    checkInfactibility();
#endif
    if ( (iout == 3) && calcFeasibilityOpt () )
      return 2;
    else if ( (iout == 1) || (iout > 2) )
      return 1;
    else
      return 0;
  }

}
