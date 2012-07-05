#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {

  /*Int Interface::vertSafeguard () {
  }*/

  // For bfgs use this
  Int Interface::vertSafeguard () {
    Vector dn (*env);
//    Real dnnorm = 0;
    Int iout = 0, ibfgs = 0;
    Bool scaleJ = ScaleVertical;

    call_ccfsg_xc (dciTrue, scaleJ);
    if (!Linear)
      this->cholesky_J ();

    NAstep (*c, dn);

//    dnnorm = dn.norm ();
    iout = dcibfgs (dn, ibfgs);
    nbfgs += ibfgs;
#ifndef NDEBUG
    checkInfactibility();
#endif
    if ( (iout == 3) && calc_feasibilityOpt () )
      return 2;
    else if ( (iout == 1) || (iout > 2) )
      return 1;
    else
      return 0;
  }

}
