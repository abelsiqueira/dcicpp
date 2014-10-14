#ifndef dci_definitions_h
#define dci_definitions_h

#include "extra.h"

#define UNUSED(x) (void)(x)

namespace DCI {

  typedef int Int;
  typedef bool Bool;
  typedef double     Real;

  typedef Int *      pInt;
  typedef Real *     pReal;
  typedef Bool *     pBool;

  const Bool dciFalse = 0;
  const Bool dciTrue = 1;
  const Real dciInf = 1e20;
  const Real dciTiny = 1e-20;
  const Real dciEps = 1e-12;

  // Unconstrained
  typedef void (*pufn)    (pInt, pInt, pReal, pReal);
  typedef void (*puofg)   (pInt, pInt, pReal, pReal, pReal, pBool);
  typedef void (*puprod)  (pInt, pInt, pBool, pReal, pReal, pReal);
  typedef void (*punames) (pInt, pInt, char*, char*);
  // Constrained
  typedef void (*pcfn)    (pInt, pInt, pInt, pReal, pReal, pReal);
  typedef void (*pcofg)   (pInt, pInt, pReal, pReal, pReal, pBool);
  typedef void (*pccfsg)  (pInt, pInt, pInt, pReal, pReal, pInt, pInt,
      pReal, pInt, pInt, pBool);
  typedef void (*pccifg)  (pInt, pInt, pInt, pReal, pReal, pReal, pBool);
  typedef void (*pcprod)  (pInt, pInt, pInt, pBool, pReal, pReal, pReal, pReal);
  typedef void (*pcnames) (pInt, pInt, pInt, char*, char*, char*);

};

#endif
