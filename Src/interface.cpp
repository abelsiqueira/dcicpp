#include "interface.h"
#include <cmath>
#include <fstream>
//#include <cassert>

namespace DCI {
  void error (int status, const char *file, int line, const char *message) {
    std::cerr << message << std::endl
              << "Status = " << status << std::endl
              << "File   = " << file << std::endl
              << "Line   = " << line << std::endl;
  }

  void Interface::assert (Bool v) {
    if (v)
      return;
    ExitFlag = -1;
    error (-1, "assert", -1, "assert failed");
    GDBSTOP ();
    throw -1;
  }

  void Interface::GDBSTOP () {
    return;
  }

  Interface::Interface () {
    //Parameters

    ufn = 0;
    cfn = 0;
    uofg = 0;
    cofg = 0;
    uprod = 0;
    cprod = 0;
    ccfsg = 0;
    ccifg = 0;
    unames = 0;
    cnames = 0;
    f = 0;
    fxc = 0;
    g = 0;
    H = 0;
    Htrip = 0;
    c = 0;
    J = 0;
    Jtrip = 0;
    LJ = 0;
    gp = 0;
    normgp = 0;
    normg = 0;
    normc = 0;
    x = 0;
    solx = 0;
    bl = 0;
    bu = 0;
    y = 0;
    yineq = 0;
    cl = 0;
    cu = 0;
    s = 0;
    sols = 0;
    xc = 0;
    sc = 0;
    feasOpt = 0;
    equatn = 0;
    linear = 0;
    nmax = 0;
    mmax = 0;
    amax = 0;
    nvar = 0;
    ncon = 0;
    nconE = 0;
    nconI = 0;
    ineqIdx = 0;
    Ineq = dciFalse;
    CurrentTime = 0;
    MaxTime = dciInf;
    DLH = dciInf;
    DLV = 0;
    Lref = dciInf;
    StartAtOne = dciFalse;
    Initialized = dciFalse;
    Running = dciFalse;
    Solved = dciFalse;
    Unbounded = dciFalse;
    UseCG = dciFalse;

//    PartialPenal = dciFalse;
    PartialPenal = dciTrue;

    project_dcp = dciTrue;
//    project_dcp = dciFalse;

//    project_dn = dciTrue;
    project_dn = dciFalse;

    project_bfgs = dciTrue;
//    project_bfgs = dciFalse;

//    trustWorstdn = dciFalse;
    trustWorstdn = dciTrue;

//    trustConvexBox = dciTrue;
    trustConvexBox = dciFalse;;

    penal_trust = dciFalse;
//    penal_trust = dciTrue;

    cholCorrection = 0;
    DisplayLevel = 1;
    env = new Environment;
    env->set_error_handler (&error);
  }

  Interface::~Interface () {
    delpointer (f);
    delpointer (fxc);
    delpointer (x);
    delpointer (solx);
    delpointer (bl);
    delpointer (bu);
    delpointer (y);
    delpointer (yineq);
    delpointer (cl);
    delpointer (cu);
    delpointer (s);
    delpointer (sols);
    delarray   (equatn);
    delarray   (ineqIdx);
    delarray   (linear);
    delpointer (xc);
    delpointer (sc);
    delpointer (feasOpt);
    delpointer (g);
    delpointer (H);
    delpointer (Htrip);
    delpointer (c);
    delpointer (gp);
    delpointer (J);
    delpointer (Jtrip);
    delpointer (LJ);
    delpointer (env); 
  }

  int Interface::start () {
    //Tests for null
    if ( (x == 0) || (bl == 0) || (bu == 0) )
      return -1;

    call_names ();
    StartTime = getTime();

    if ( (ncon > 0) && (equatn == 0) ) {
      equatn = new Bool[ncon];
      for (Int i = 0; i < ncon; i++)
        equatn[i] = dciTrue;
      nconE = ncon;
      nconI = 0;
      Ineq = dciFalse;
    }

    if (linear == 0) {
      linear = new Bool[ncon];
      for (Int i = 0; i < ncon; i++)
        linear[i] = dciFalse;
      nconNL = ncon;
      nconL = 0;
      Linear = dciFalse;
    }

    if (nmax == 0)
      nmax = nvar;
    if (mmax == 0)
      mmax = ncon;
    if (amax == 0) {
      amax = nvar*ncon;
      if (amax < 10)
        amax = 10;
    }

    amax = amax + nconI;

    f = new Real;
    fxc = new Real;
    g = new Vector (*env, nmax + nconI);
    gx = g->get_doublex();
    if (ncon > 0) {
      c = new Vector (*env, mmax);
      cx = c->get_doublex();
      if (Ineq) {
        s = new Vector (*env, nconI);
        sx = s->get_doublex();
        sc = new Vector (*env, nconI);
        scx = sc->get_doublex();
      }
      Jtrip = new Triplet (*env, mmax, nmax + nconI, amax, 0); //0 is for unsymmetric
      Jx = Jtrip->get_doublex();
      Ji = static_cast < Int * > (Jtrip->triplet->i);
      Jj = static_cast < Int * > (Jtrip->triplet->j);
      for (Int i = 0; i < amax; i++) {
        Jx[i] = 0;
        Ji[i] = 0;
        Jj[i] = 0;
      }
      LJ = new Factor (*env);
    }
    yineq = new Vector (*env, nvar + nconI);
    yineqx = yineq->get_doublex ();
    gp = new Vector (*env, nmax + nconI);


#ifdef LOCALTEST
    if (solx != 0) {
      Vector tmp (*xc);
      *xc = *solx;
      call_fn_xc ();
      sols = new Vector (*env, nconI);
      pReal ps = sols->get_doublex ();
      for (Int i = 0; i < nconI; i++) {
        ps[i] = cx[ineqIdx[i]];
      }
    }
#endif

    InitialParameters ();
    InitialValues ();
    Initialized = dciTrue;

    for (Int i = 0; i < nvar; i++) {
      assert (blx[i] < xx[i]);
      assert (xx[i] < bux[i]);
    }

    return 0;
  }

  void Interface::show (Int displvl) {
    show (std::cout, displvl);
  }

  void Interface::show (std::ostream & out, Int displvl) {
    if (!Initialized) {
      out << "Problem not initialized" << std::endl;
      return;
    }

    Real yoff = 0;
    for (Int i = 0; i < nconI; i++) {
      Real yi = yx[ineqIdx[i]], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      if (li > -dciInf)
        yoff += Max(yi, 0.0);
      if (ui < dciInf)
        yoff += Max(-yi, 0.0);
    }

    DisplayLevel = displvl;

    out << "**********************************************************" << std::endl;
    out << "Problem name: " << problemName << std::endl
        << std::endl
        << "Number of Variables: " << nvar << std::endl
        << "Number of Constraints: " << ncon << std::endl 
        << "             Equality: " << nconE << std::endl
        << "           Inequality: " << nconI << std::endl
        << std::endl;

    if (!Solved) {
      out << "Problem not solved yet" << std::endl;
    } else {
      if (ExitFlag == 0)
        out << "The Algorithm has Converged" << std::endl;
      else if (ExitFlag == 1)
        out << "rhomax became too short" << std::endl;
      else if (ExitFlag == 2)
        out << "The maximum number of iterations was reached" << std::endl;
      else if (ExitFlag == 3)
        out << "The maximum number of restorarions was reached" << std::endl;
      else if (ExitFlag == 4)
        out << "The restoration has failed" << std::endl;
      else if (ExitFlag == 5)
        out << "The step became too short" << std::endl;
      else if (ExitFlag == 6)
        out << "The problem is unbounded" << std::endl;
      else if (ExitFlag == 7)
        out << "The problem reached the time limit" << std::endl;
      else if (ExitFlag == 8)
        out << "The problem infeasible" << std::endl;

      out << "f(x) = " << *f << std::endl
          << "|c(x)| = " << normc << std::endl
          << "|g(x) + J(x)'*y| = " << normgp << std::endl
          << "y offset = " << yoff << std::endl
          << "BFGS? " << ((tbfgs > 0) ? "yes" : "no") << std::endl
          << "Number of Iterations = " << iter << std::endl
          << "Elapsed Time = " << (CurrentTime > 0 ? CurrentTime : 0) << " s" << std::endl;

      if (DisplayLevel > 2) {
        std::cout << std::endl
            << "----------------------" << std::endl
            << "  Iter   = " << iter << std::endl
            << "  f(x)   = " << *f << std::endl
            << "  |c(x)| = " << c->norm () << std::endl
            << "  |gp|   = " << gp->norm () << std::endl
            << "  |g|    = " << g->norm () << std::endl
            << "  ngp    = " << ngp << std::endl
            << "  rho    = " << rho << std::endl
            << "  rhomax = " << rhomax << std::endl
            << "  mu     = " << mu << std::endl
            << "  DLH    = " << DLH << std::endl
            << "  DLV    = " << DLV << std::endl
            << "  DeltaH = " << DeltaH << std::endl
            << "  DeltaV = " << DeltaV << std::endl
            << "  tSoc   = " << tSoc << std::endl
            << "  tSteih = " << tSteih << std::endl
            << "  tRej   = " << tRej << std::endl
            << "  tbfgs  = " << tbfgs << std::endl
            << "  tRest  = " << tRest << std::endl
            << "  itssmll= " << itssmll << std::endl
            << "  gap  = " << gap << std::endl
            << "----------------------" << std::endl;
        if (ncon > 0) {
          std::cout << "lambda = " << std::endl;
          y->print_more();
          if (Ineq) {
            std::cout << "s = " << std::endl;
            s->print_more();
            std::cout << "sc = " << std::endl;
            sc->print_more();
          }
        }
      }
      if (DisplayLevel > 1) {
        if (nvar <= 10) {
          out << "x: " << std::endl;
          x->print_more ();
        } else
          out << "Number of Variables is too big to show: Increase nvarshowmax" << std::endl;
        if ( (ncon > 0) && (ncon <= 10) ) {
          out << "c: " << std::endl;
          c->print_more ();
          out << "y: " << std::endl;
          y->print_more ();
        } else if (ncon > 0)
          out << "Number of Constraints is too big to show: Increase nconshowmax" << std::endl;
      }
    }
    out << "**********************************************************" << std::endl;

  }

  // Problem name & nvar & ncon & iters & time & converged
  void Interface::printLatex (char * filename) const {
    std::ofstream file;
    if (filename == 0) {
      switch (ExitFlag) {
        case -1:
          file.open ("latex_assert", std::ios_base::app);
          break;
        case 0:
          file.open ("latex_convergence", std::ios_base::app);
          break;
        case 1:
          file.open ("latex_rhomax", std::ios_base::app);
          break;
        case 2:
          file.open ("latex_maxiter", std::ios_base::app);
          break;
        case 3:
          file.open ("latex_maxrest", std::ios_base::app);
          break;
        case 4:
          file.open ("latex_restfail", std::ios_base::app);
          break;
        case 5:
          file.open ("latex_shortstep", std::ios_base::app);
          break;
        case 6:
          file.open ("latex_unbounded", std::ios_base::app);
          break;
        case 7:
          file.open ("latex_timelimit", std::ios_base::app);
          break;
      }
    } else
      file.open (filename, std::ios_base::app);
    file << problemName << " & "
         << nvar << " & "
         << ncon << " & "
         << iter << " & "
         << (CurrentTime > 0 ? CurrentTime : 0) << " & "
         << ((ncon > 0) ? "con" : "unc") << " & "
         << ((tbfgs > 0) ? "bfgs" : "") << "\\\\ \\hline\n";

    file.close ();
  }

  void Interface::set_x (size_t n, Real * V) {
    x = new Vector (*env, n, V);
    set_nvar(n);
    xx = x->get_doublex();
    xc = new Vector (*env, n, V);
    xcx = xc->get_doublex();
  }

  void Interface::set_sol (size_t n, Real * V) {
    solx = new Vector (*env, n, V);
  }

  void Interface::set_bl (size_t n, Real * V) {
    bl = new Vector (*env, n, V);
    set_nvar(n);
    blx = bl->get_doublex();
  }

  void Interface::set_bu (size_t n, Real * V) {
    bu = new Vector (*env, n, V);
    set_nvar(n);
    bux = bu->get_doublex();
  }

  void Interface::set_lambda (size_t n, Real * V) {
    y = new Vector (*env, n, V);
    set_ncon(n);
    yx = y->get_doublex();
  }

  void Interface::set_cl (size_t n, Real * V) {
    cl = new Vector (*env, n, V);
    set_ncon(n);
    clx = cl->get_doublex();
  }

  void Interface::set_cu (size_t n, Real * V) {
    cu = new Vector (*env, n, V);
    set_ncon(n);
    cux = cu->get_doublex();
  }

  void Interface::set_equatn (size_t n, Bool * V) {
    equatn = new Bool[n];
    set_ncon(n);
    nconE = nconI = 0;
    for (size_t i = 0; i < n; i++) {
      Bool tmp = V[i];
      equatn[i] = tmp;
      if (tmp == dciFalse)
        nconI++;
      else
        nconE++;
    }
    if (nconI > 0) {
      Ineq = dciTrue;
      ineqIdx = new Int[nconI];
      Int numI = 0;
      for (size_t i = 0; i < n; i++) {
        if (equatn[i] == dciFalse) {
          ineqIdx[numI] = i;
          numI++;
        }
      }
    }
  }

  void Interface::set_linear (size_t n, Bool * V) {
    linear = new Bool[n];
    set_ncon(n);
    nconL = nconNL = 0;
    for (size_t i = 0; i < n; i++) {
      Bool tmp = V[i];
      linear[i] = tmp;
      if (tmp == dciFalse)
        nconNL++;
      else
        nconL++;
    }
    if (nconNL == 0)
      Linear = dciTrue;
  }

  void Interface::call_fn () {
    if (ncon == 0) {
      (*ufn) (&nvar, xx, f);
    } else {
      (*cfn) (&nvar, &ncon, xx, f, &mmax, cx);
    }
    if (Running) {
      *f -= mu*calc_pen ();
      Int numI = 0;
      for (Int i = 0; i < ncon; i++) {
        if (equatn[i] == dciFalse) {
          cx[i] -= sx[numI];
          numI++;
        }
      }
    }
  }

  void Interface::call_fn_xc () {
    if (ncon == 0) {
      (*ufn) (&nvar, xcx, f);
    } else {
      (*cfn) (&nvar, &ncon, xcx, fxc, &mmax, cx);
    }
    if (Running) {
      *fxc -= mu*calc_pen_xc ();
      Int numI = 0;
      for (Int i = 0; i < ncon; i++) {
        if (equatn[i] == dciFalse) {
          cx[i] -= scx[numI];
          numI++;
        }
      }
    }
  }

  void Interface::call_ofg (Bool grad) {
    if (ncon == 0)
      (*uofg) (&nvar, xx, f, gx, &grad);
    else
      (*cofg) (&nvar, xx, f, gx, &grad);
    if (Running) {
      Real val = 0.0;
      for (Int i = 0; i < nvar; i++) {
        Real xi = xx[i], bli = blx[i], bui = bux[i];
        if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (PartialPenal) {
            if ( (xi - bli) < (bui - xi) ) {
              gx[i] *= (xi - bli);
              gx[i] -= mu;
              val += log (xi - bli);
            } else {
              gx[i] *= (bui - xi);
              gx[i] += mu;
              val += log (bui - xi);
            }
          } else {
            gx[i] *= (xi - bli) * (bui - xi);
            gx[i] += mu * (2*xi - bli - bui);
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          gx[i] *= (bui - xi);
          gx[i] += mu;
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          gx[i] *= (xi - bli);
          gx[i] -= mu;
          val += log (xi - bli);
        }
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = sx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (PartialPenal) {
            if ( (si - cli) < (cui - si) ) {
              gx[j] = -mu;
              val += log (si - cli);
            } else {
              gx[j] = mu;
              val += log (cui - si);
            }
          } else {
            gx[j] = mu * (2*si - cli - cui);
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          gx[j] = mu;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          gx[j] = -mu;
          val += log (si - cli);
        }
      }
      *f -= mu*val;
    }
  }

  void Interface::call_ofg_xc (Bool grad) {
    if (ncon == 0)
      (*uofg) (&nvar, xcx, fxc, gx, &grad);
    else
      (*cofg) (&nvar, xcx, fxc, gx, &grad);
    if (Running) {
      Real val = 0.0;
      for (Int i = 0; i < nvar; i++) {
        Real xi = xcx[i], bli = blx[i], bui = bux[i];
        if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (PartialPenal) {
            if ( (xi - bli) < (bui - xi) ) {
              gx[i] *= (xi - bli);
              gx[i] -= mu;
              val += log (xi - bli);
            } else {
              gx[i] *= (bui - xi);
              gx[i] += mu;
              val += log (bui - xi);
            }
          } else {
            gx[i] *= (xi - bli) * (bui - xi);
            gx[i] += mu * (2*xi - bli - bui);
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          gx[i] *= (bui - xi);
          gx[i] += mu;
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          gx[i] *= (xi - bli);
          gx[i] -= mu;
          val += log (xi - bli);
        }
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = scx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (PartialPenal) {
            if ( (si - cli) < (cui - si) ) {
              gx[j] = -mu;
              val += log (si - cli);
            } else {
              gx[j] = mu;
              val += log (cui - si);
            }
          } else {
            gx[j] = mu * (2*si - cli - cui);
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          gx[j] = mu;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          gx[j] = -mu;
          val += log (si - cli);
        }
      }
      *fxc -= mu*val;
    }
  }

  void Interface::call_ccfsg (Bool grad, Bool scale) {
    pInt nnzj = new Int(0);
    (*ccfsg) (&nvar, &ncon, xx, &mmax, cx, nnzj, &amax, Jx, Jj, Ji, &grad);
    if (grad == dciTrue) {
      if (StartAtOne) {
        for (Int i = 0; i < *nnzj; i++) {
          Jj[i]--;
          Ji[i]--;
        }
      }

      if (Running) {
        Vector Diag(*env);
        Diag.reset (nvar + nconI, 1);
        pReal Diagx = Diag.get_doublex();
        if (scale) { 
          scale_x (Diag);
          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jj[k];
            Jx[k] *= Diagx[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = sx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -Diagx[j];
            else
              Jx[*nnzj] = -1;
            Jj[*nnzj] = j;
            Ji[*nnzj] = i;
            (*nnzj)++;
            numI++;
            j++;
          }
        }
      }

      *(Jtrip->get_pnnz()) = *nnzj;
      delpointer (J);
      J = new Sparse (*env);
      J->triplet_to_sparse (*Jtrip, amax);
    } else {
      if (Running) {
        Int numI = 0;
        for (Int i = 0; i < ncon; i++) {
          if (equatn[i] == dciFalse) {
            Real si = sx[numI];
            cx[i] -= si;
            numI++;
          }
        }
      }
    }
    delete nnzj;
  }
  
  void Interface::call_ccfsg_xc (Bool grad, Bool scale) {
    pInt nnzj = new Int(0);
    (*ccfsg) (&nvar, &ncon, xcx, &mmax, cx, nnzj, &amax, Jx, Jj, Ji, &grad);
    if (grad == dciTrue) {
      if (StartAtOne) {
        for (Int i = 0; i < *nnzj; i++) {
          Jj[i]--;
          Ji[i]--;
        }
      }

      if (Running) {
        Vector Diag(*env);
        Diag.reset (nvar + nconI, 1);
        pReal Diagx = Diag.get_doublex();
        if (scale) {
          scale_xc (Diag);

          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jj[k];
            Jx[k] *= Diagx[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = scx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -Diagx[j];
            else
              Jx[*nnzj] = -1;
            Jj[*nnzj] = j;
            Ji[*nnzj] = i;
            (*nnzj)++;
            numI++;
            j++;
          }
        }
      }

      *(Jtrip->get_pnnz()) = *nnzj;
      delpointer (J);
      J = new Sparse (*env);
      J->triplet_to_sparse (*Jtrip, amax);
    } else {
      if (Running) {
        Int numI = 0;
        for (Int i = 0; i < ncon; i++) {
          if (equatn[i] == dciFalse) {
            Real si = scx[numI];
            cx[i] -= si;
            numI++;
          }
        }
      }
    }
    delete nnzj;
  }

  void Interface::call_prod (Bool gotder, pReal px, pReal ux) {
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      Real xi = xx[i], bli = blx[i], bui = bux[i], pxi = px[i];
      if ( (bli > -dciInf) && (bui < dciInf) ) {
        if (PartialPenal) {
          if ( (xi - bli) < (bui - xi) )
            ppx[i] = (xi - bli) * pxi;
          else
            ppx[i] = (bui - xi) * pxi;
        } else
          ppx[i] = (xi - bli) * (bui - xi) * pxi;
      } else if ( (bli <= -dciInf) && (bui < dciInf) )
        ppx[i] = (bui - xi) * pxi;
      else if ( (bli > -dciInf) && (bui >= dciInf) )
        ppx[i] = (xi - bli) * pxi;
      else
        ppx[i] = pxi;
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xx, &mmax, yx, ppx, ux);
    for (Int i = 0; i < nvar; i++) {
      Real xi = xx[i], bli = blx[i], bui = bux[i], pxi = px[i], uxi = ux[i];
      if ( (bli > -dciInf) && (bui < dciInf) ) {
        if (PartialPenal) {
          if ( (xi - bli) < (bui - xi) )
            ux[i] = (xi - bli) * uxi + mu * pxi;
          else
            ux[i] = (bui - xi) * uxi + mu * pxi;
        } else
          ux[i] = (xi - bli) * (bui - xi) * uxi + mu * ( (xi - bli)*(xi - bli) + (xi - bui)*(xi - bui) ) * pxi;
      } else if ( (bli <= -dciInf) && (bui < dciInf) )
        ux[i] = (bui - xi) * uxi + mu * pxi;
      else if ( (bli > -dciInf) && (bui >= dciInf) )
        ux[i] = (xi - bli) * uxi + mu * pxi;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real si = sx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], pxi = px[j];
      if ( (cli > -dciInf) && (cui < dciInf) ) {
        if (PartialPenal) {
          if ( (si - cli) < (cui - si) )
            ux[j] = Max(mu, minBk) * pxi;
          else
            ux[j] = Max(mu, minBk) * pxi;
        } else
          ux[j] = Max(mu,minBk) * ( (si - cli)*(si - cli) + (si - cui)*(si - cui) ) * pxi;
      } else if ( (cli <= -dciInf) && (cui < dciInf) )
        ux[j] = Max(mu,minBk) * pxi;
      else if ( (cli > -dciInf) && (cui >= dciInf) )
        ux[j] = Max(mu,minBk) * pxi;
      else
        ux[j] = minBk;
    }
  }

  void Interface::call_prod_xc (Bool gotder, pReal px, pReal ux) {
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], bli = blx[i], bui = bux[i], pxi = px[i];
      if ( (bli > -dciInf) && (bui < dciInf) ) {
        if (PartialPenal) {
          if ( (xi - bli) < (bui - xi) )
            ppx[i] = (xi - bli) * pxi;
          else
            ppx[i] = (bui - xi) * pxi;
        } else
          ppx[i] = (xi - bli) * (bui - xi) * pxi;
      } else if ( (bli <= -dciInf) && (bui < dciInf) )
        ppx[i] = (bui - xi) * pxi;
      else if ( (bli > -dciInf) && (bui >= dciInf) )
        ppx[i] = (xi - bli) * pxi;
      else
        ppx[i] = pxi;
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xcx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xcx, &mmax, yx, ppx, ux);
    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], bli = blx[i], bui = bux[i], pxi = px[i], uxi = ux[i];
      if ( (bli > -dciInf) && (bui < dciInf) ) {
        if (PartialPenal) {
          if ( (xi - bli) < (bui - xi) )
            ux[i] = (xi - bli) * uxi + mu * pxi;
          else
            ux[i] = (bui - xi) * uxi + mu * pxi;
        } else
          ux[i] = (xi - bli) * (bui - xi) * uxi + mu * ( (xi - bli)*(xi - bli) + (xi - bui)*(xi - bui) ) * pxi;
      } else if ( (bli <= -dciInf) && (bui < dciInf) )
        ux[i] = (bui - xi) * uxi + mu * pxi;
      else if ( (bli > -dciInf) && (bui >= dciInf) )
        ux[i] = (xi - bli) * uxi + mu * pxi;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real si = scx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], pxi = px[j];
      if ( (cli > -dciInf) && (cui < dciInf) ) {
        if (PartialPenal)
          ux[j] = Max(mu, minBk) * pxi;
        else
          ux[j] = Max(mu,minBk) * ( (si - cli)*(si - cli) + (si - cui)*(si - cui) ) * pxi;
      } else if ( (cli <= -dciInf) && (cui < dciInf) )
        ux[j] = Max(mu,minBk) * pxi;
      else if ( (cli > -dciInf) && (cui >= dciInf) )
        ux[j] = Max(mu,minBk) * pxi;
      else
        ux[j] = minBk;
    }
  }

  void Interface::call_names () {
    if ( (unames == 0) && (cnames == 0) ) {
      problemName = "No name";
      return;
    }
    char pname[10], vnames[10*nvar];
    if (ncon == 0)
      (*unames) (&nvar, pname, vnames);
    else {
      char gnames[10*ncon];
      (*cnames) (&nvar, &ncon, pname, vnames, gnames);
    }
    problemName = pname;
  }

  void Interface::analyze_J () {
    if (LJ == 0) 
      LJ = new Factor (*env);
    LJ->analyze (*J);
  }

  void Interface::cholesky_J () {
    if (LJ == 0)
      std::cerr << "analyze should be called first" << std::endl;
//    if (cholCorrection > 0)
//      cholCorrection = normc;
    LJ->factorize (*J, cholCorrection);
    cholFacs++;
    if (!env->IsPosDef()) {
      cholCorrection = 1;
      LJ->factorize (*J, cholCorrection);
      cholFacs++;
    }
  }

  Real Interface::getTime () {
    clock_t t = clock();

    return static_cast<Real>(static_cast<Real>(t)/CLOCKS_PER_SEC);
  }

  void Interface::checkInfactibility () {
    for (Int i = 0; i < nvar; i++) {
      if ( (xx[i] > 1e10) || (xx[i] < -1e10) || (xcx[i] > 1e10) || (xcx[i] < -1e10) ) {
        Unbounded = dciTrue;
        break;
      }
      if (xx[i] >= bux[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xx[i] < bux[i]);

      if (xx[i] <= blx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xx[i] > blx[i]);

      if (xcx[i] >= bux[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xcx[i] < bux[i]);

      if (xcx[i] <= blx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xcx[i] > blx[i]);
    }
    if (Unbounded)
      return;
    for (Int i = 0; i < nconI; i++) {
      if ( (sx[ineqIdx[i]] > 1e10) || (sx[ineqIdx[i]] < -1e10) || (scx[ineqIdx[i]] > 1e10) || (scx[ineqIdx[i]] < -1e10) ) {
        Unbounded = dciTrue;
        break;
      }
      assert (sx[i] < cux[ineqIdx[i]]);
      assert (sx[i] > clx[ineqIdx[i]]);
      assert (scx[i] < cux[ineqIdx[i]]);
      assert (scx[i] > clx[ineqIdx[i]]);
    }
  }


}
