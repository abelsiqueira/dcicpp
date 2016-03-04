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
  typedef void (*pufn)    (pInt, const pInt, const pReal, pReal);
  typedef void (*puofg)   (pInt, const pInt, const pReal, pReal, pReal, const pBool);
  typedef void (*puprod)  (pInt, const pInt, const pBool, const pReal, const pReal, pReal);
  typedef void (*punames) (pInt, const pInt, char*, char*);
  // Constrained
  typedef void (*pcfn)    (pInt, const pInt, const pInt, const pReal, pReal, pReal);
  typedef void (*pcofg)   (pInt, const pInt, const pReal, pReal, pReal, pBool);
  typedef void (*pccfsg)  (pInt, const pInt, const pInt, const pReal, pReal, pInt, const pInt,
      pReal, pInt, pInt, const pBool);
  typedef void (*pccifg)  (pInt, const pInt, const pInt, const pReal, pReal, pReal, const pBool);
  typedef void (*pcprod)  (pInt, const pInt, const pInt, const pBool, const pReal, const pReal, pReal, pReal);
  typedef void (*pcnames) (pInt, const pInt, const pInt, char*, char*, char*);

};

#endif
