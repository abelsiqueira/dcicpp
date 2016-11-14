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
  typedef void (*pufn)    (pInt, const Int *, const Real *, pReal);
  typedef void (*puofg)   (pInt, const Int *, const Real *, pReal, pReal, const Bool *);
  typedef void (*puprod)  (pInt, const Int *, const Bool *, const Real *, const Real *, pReal);
  typedef void (*punames) (pInt, const Int *, char*, char*);
  // Constrained
  typedef void (*pcfn)    (pInt, const Int *, const Int *, const Real *, pReal, pReal);
  typedef void (*pcofg)   (pInt, const Int *, const Real *, pReal, pReal, pBool);
  typedef void (*pccfsg)  (pInt, const Int *, const Int *, const Real *, pReal, pInt, const Int *,
      pReal, pInt, pInt, const Bool *);
  typedef void (*pccifg)  (pInt, const Int *, const Int *, const Real *, pReal, pReal, const Bool *);
  typedef void (*pcprod)  (pInt, const Int *, const Int *, const Bool *, const Real *, const Real *, pReal, pReal);
  typedef void (*pcnames) (pInt, const Int *, const Int *, char*, char*, char*);

};

#endif
