#ifndef dci_definitions_h
#define dci_definitions_h

#include "extra.h"

namespace DCI {

#ifdef LONG
  typedef long int   Int;
  typedef long int   Bool;
#else
  typedef int   Int;
  typedef long int   Bool;
#endif

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
  typedef void (*pufn)    (pInt, pReal, pReal);
  typedef void (*puofg)   (pInt, pReal, pReal, pReal, pBool);
  typedef void (*puprod)  (pInt, pBool, pReal, pReal, pReal);
  typedef void (*punames) (pInt, char*, char*);
  // Constrained
  typedef void (*pcfn)    (pInt, pInt, pReal, pReal, pInt, pReal);
  typedef void (*pcofg)   (pInt, pReal, pReal, pReal, pBool);
  typedef void (*pccfsg)  (pInt, pInt, pReal, pInt, pReal, pInt, pInt, pReal, pInt, pInt, pBool);
  typedef void (*pccifg)  (pInt, pInt, pReal, pReal, pReal, pBool);
  typedef void (*pcprod)  (pInt, pInt, pBool, pReal, pInt, pReal, pReal, pReal);
  typedef void (*pcnames) (pInt, pInt, char*, char*, char*);

};

#endif
