#include "interface.h"
#include <cmath>
#include <sstream>
#include <algorithm>
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
    throw "assert failed";
  }

  void Interface::GDBSTOP () {
    return;
  }

  Interface::Interface () {
    //Parameters
    Initialization ();
    //Mumps
    if (UseMUMPS) {
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
      id.job=JOB_INIT;
      id.par=1;
      id.sym=2;
      id.comm_fortran=USE_COMM_WORLD;
      dmumps_c(&id);
    }

//    cholCorrection = 1e-4;
    cholCorrection = 0;
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
    delpointer (Lambda);
    delpointer (initial_x);
  }

  int Interface::start () {
    //Tests for null
    if ( (x == 0) || (bl == 0) || (bu == 0) )
      return -1;

    //Mumps
    id.n = ncon + nvar + nconI;

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
    if (UseMUMPS)
      amax += nvar + nconI + ncon;

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

    DefineParameters ();
    InitialValues ();
    Initialized = dciTrue;

    for (Int i = 0; i < nvar; i++) {
      assert (blx[i] < xx[i]);
      assert (xx[i] < bux[i]);
    }

    return 0;
  }

  void Interface::show (std::ostream & out) {
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

    if (DisplayLevel > 0) {
      out << "**********************************************************" << std::endl;
      out << "Problem name: " << problemName << std::endl
          << std::endl
          << "Number of Variables: " << nvar << std::endl
          << "Number of Constraints: " << ncon << std::endl 
          << "             Equality: " << nconE << std::endl
          << "           Inequality: " << nconI << std::endl
          << std::endl;
    }

#ifdef PRINT_MATLAB
    if (J != 0) {
      std::string matlab_filename("matrix_jacob");
      matlab_filename += problemName;
      switch (ExitFlag) {
        case -1:
          matlab_filename += "_assert";
          break;
        case 0:
          matlab_filename += "_conv";
          break;
        case 1:
          matlab_filename += "_rhomax";
          break;
        case 2:
          matlab_filename += "_maxiter";
          break;
        case 3:
          matlab_filename += "_maxrest";
          break;
        case 4:
          matlab_filename += "_restfail";
          break;
        case 5:
          matlab_filename += "_short";
          break;
        case 6:
          matlab_filename += "_unlimited";
          break;
        case 7:
          matlab_filename += "_time";
          break;
        case 8:
          matlab_filename += "_infeasible";
          break;
      }
      matlab_filename += ".m";
      matlab_filename.erase( remove_if(matlab_filename.begin(), 
            matlab_filename.end(), isspace), matlab_filename.end());
      std::ofstream matlab_file(matlab_filename.c_str());
      J->print_matlab(matlab_file);
      matlab_file.close();
    }
#endif

    if (!Solved) {
      out << "Problem not solved yet" << std::endl;
    } else {
      out << "EXIT: ";
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
        out << "The problem is unlimited" << std::endl;
      else if (ExitFlag == 7)
        out << "The problem reached the time limit" << std::endl;
      else if (ExitFlag == 8)
        out << "Stopped at a stationary for the infeasibility" << std::endl;

      if (DisplayLevel > 0) {
        out << "f(x) = " << *f << std::endl
            << "|c(x)| = " << normc << std::endl
            << "|g(x) + J(x)'*y| = " << normgp << std::endl
            << "y offset = " << yoff << std::endl
            << "BFGS? " << ((tbfgs > 0) ? "yes" : "no") << std::endl
            << "Number of Iterations = " << iter << std::endl
            << "Elapsed Time = " << (CurrentTime > 0 ? CurrentTime : 0) << " s" << std::endl;
      }

      if (DisplayLevel > 2) {
        std::cout << std::endl
            << "----------------------" << std::endl
            << "  Iter   = " << iter << std::endl
            << "  f(x)   = " << *f << std::endl
            << "  |c(x)| = " << (ncon > 0 ? c->norm () : 0) << std::endl
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
          file.open ("latex_unlimited", std::ios_base::app);
          break;
        case 7:
          file.open ("latex_timelimit", std::ios_base::app);
          break;
        case 8:
          file.open ("latex_infeasible", std::ios_base::app);
          break;
        default:
          std::stringstream aux;
          aux << "Exitflag value " << ExitFlag;
          throw(aux.str());
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
    initial_x = new Real[n];
    for (Int i = 0; i < (Int)n; i++)
      initial_x[i] = 1.0;
  }

  void Interface::set_sol (size_t n, Real * V) {
    solx = new Vector (*env, n, V);
  }

  void Interface::set_bl (size_t n, Real * V) {
    bl = new Vector (*env, n, V);
    set_nvar(n);
    blx = bl->get_doublex();
    for (size_t i = 0; i < n; i++) {
      if (blx[i] > -dciInf) {
        Bounded = dciTrue;
        break;
      }
    }
  }

  void Interface::set_bu (size_t n, Real * V) {
    bu = new Vector (*env, n, V);
    set_nvar(n);
    bux = bu->get_doublex();
    for (size_t i = 0; i < n; i++) {
      if (bux[i] < dciInf) {
        Bounded = dciTrue;
        break;
      }
    }
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
      Bounded = dciTrue;
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
      *f /= objfun_scale;
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
      *fxc /= objfun_scale;
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
      if (objfun_scale != 1) {
        *f /= objfun_scale;
        if (grad) {
          for (Int i = 0; i < nvar; i++)
            gx[i] /= objfun_scale;
        }
      }
      for (Int i = 0; i < nvar; i++) {
        Real xi = xx[i], bli = blx[i], bui = bux[i];
        if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (PartialPenal) {
            if ( (xi - bli) < (bui - xi) ) {
              if (grad) {
                gx[i] *= (xi - bli);
                gx[i] -= mu;
              }
              val += log (xi - bli);
            } else {
              if (grad) {
                gx[i] *= (bui - xi);
                gx[i] += mu;
              }
              val += log (bui - xi);
            }
          } else {
            if (grad) {
              gx[i] *= (xi - bli) * (bui - xi);
              gx[i] += mu * (2*xi - bli - bui);
            }
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          if (grad) {
            gx[i] *= (bui - xi);
            gx[i] += mu;
          }
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          if (grad) {
            gx[i] *= (xi - bli);
            gx[i] -= mu;
          }
          val += log (xi - bli);
        }
        if (grad)
          gx[i] *= initial_x[i];
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = sx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (PartialPenal) {
            if ( (si - cli) < (cui - si) ) {
              if (grad)
                gx[j] = -mu;
              val += log (si - cli);
            } else {
              if (grad)
                gx[j] = mu;
              val += log (cui - si);
            }
          } else {
            if (grad)
              gx[j] = mu * (2*si - cli - cui);
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          if (grad)
            gx[j] = mu;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          if (grad)
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
      if (objfun_scale != 1) {
        *fxc /= objfun_scale;
        if (grad) {
          for (Int i = 0; i < nvar; i++)
            gx[i] /= objfun_scale;
        }
      }
      for (Int i = 0; i < nvar; i++) {
        Real xi = xcx[i], bli = blx[i], bui = bux[i];
        if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (PartialPenal) {
            if ( (xi - bli) < (bui - xi) ) {
              if (grad) {
                gx[i] *= (xi - bli);
                gx[i] -= mu;
              }
              val += log (xi - bli);
            } else {
              if (grad) {
                gx[i] *= (bui - xi);
                gx[i] += mu;
              }
              val += log (bui - xi);
            }
          } else {
            if (grad) {
              gx[i] *= (xi - bli) * (bui - xi);
              gx[i] += mu * (2*xi - bli - bui);
            }
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          if (grad) {
            gx[i] *= (bui - xi);
            gx[i] += mu;
          }
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          if (grad) {
            gx[i] *= (xi - bli);
            gx[i] -= mu;
          }
          val += log (xi - bli);
        }
        if (grad)
          gx[i] *= initial_x[i];
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = scx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (PartialPenal) {
            if ( (si - cli) < (cui - si) ) {
              if (grad)
                gx[j] = -mu;
              val += log (si - cli);
            } else {
              if (grad)
                gx[j] = mu;
              val += log (cui - si);
            }
          } else {
            if (grad)
              gx[j] = mu * (2*si - cli - cui);
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          if (grad)
            gx[j] = mu;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          if (grad)
            gx[j] = -mu;
          val += log (si - cli);
        }
      }
      *fxc -= mu*val;
    }
  }

  void Interface::call_ccfsg (Bool grad, Bool scale) {
    UpdateScaling_x();
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
        if (scale) { 
          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jj[k];
            Jx[k] *= Lambda[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = sx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -Lambda[j];
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
      //Mumps
      if (UseMUMPS) {
        for (Int i = 0; i < *nnzj; i++) {
          Ji[i] += nvar + nconI + 1;
          Jj[i]++;
        }
        for (Int i = 0; i < nvar + nconI; i++) {
          Ji[*nnzj + i] = i + 1;
          Jj[*nnzj + i] = i + 1;
          Jx[*nnzj + i] = 1;
        }
        if (cholCorrection > 0) {
          for (Int i = 0; i < ncon; i++) {
            Ji[*nnzj + nvar + nconI + i] = nvar + nconI + i + 1;
            Jj[*nnzj + nvar + nconI + i] = nvar + nconI + i + 1;
            Jx[*nnzj + nvar + nconI + i] = cholCorrection;
          }
        }
        id.irn = reinterpret_cast<int*>(Ji);
        id.jcn = reinterpret_cast<int*>(Jj);
        id.a = reinterpret_cast<double*>(Jx);
        if (cholCorrection > 0)
          id.nz = (int)(*nnzj) + nvar + nconI + ncon;
        else
          id.nz = (int)(*nnzj) + nvar + nconI;
      }
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
    UpdateScaling_xc();
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
        if (scale) {
          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jj[k];
            Jx[k] *= Lambda[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = scx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -Lambda[j];
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
      //Mumps
      if (UseMUMPS) {
        for (Int i = 0; i < *nnzj; i++) {
          Ji[i] += nvar + nconI + 1;
          Jj[i]++;
        }
        for (Int i = 0; i < nvar + nconI; i++) {
          Ji[*nnzj + i] = i + 1;
          Jj[*nnzj + i] = i + 1;
          Jx[*nnzj + i] = 1;
        }
        if (cholCorrection > 0) {
          for (Int i = 0; i < ncon; i++) {
            Ji[*nnzj + nvar + nconI + i] = nvar + nconI + i + 1;
            Jj[*nnzj + nvar + nconI + i] = nvar + nconI + i + 1;
            Jx[*nnzj + nvar + nconI + i] = cholCorrection;
          }
        }
        id.irn = reinterpret_cast<int*>(Ji);
        id.jcn = reinterpret_cast<int*>(Jj);
        id.a = reinterpret_cast<double*>(Jx);
        if (cholCorrection > 0)
          id.nz = (int)(*nnzj) + nvar + nconI + ncon;
        else
          id.nz = (int)(*nnzj) + nvar + nconI;
      }
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
    UpdateScaling_x();
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      ppx[i] = Lambda[i] * px[i];
    }
    if (objfun_scale != 1) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= objfun_scale;
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xx, &mmax, yx, ppx, ux);
    if (objfun_scale != 1) {
      for (Int i = 0; i < nvar; i++)
        ux[i] /= objfun_scale;
      for (Int i = 0; i < ncon; i++)
        yx[i] /= objfun_scale;
    }
    for (Int i = 0; i < nvar; i++) {
      ux[i] *= Lambda[i];
      Real li = blx[i], ui = bux[i], pxi = px[i];
      if (li > -dciInf || ui < dciInf)
        ux[i] += mu * pxi;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      ux[j] = Max(mu, minBk) * px[j];
    }
  }

  void Interface::call_prod_xc (Bool gotder, pReal px, pReal ux) {
    UpdateScaling_xc();
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      ppx[i] = Lambda[i] * px[i];
    }
    if (objfun_scale != 1) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= objfun_scale;
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xcx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xcx, &mmax, yx, ppx, ux);
    if (objfun_scale != 1) {
      for (Int i = 0; i < nvar; i++)
        ux[i] /= objfun_scale;
      for (Int i = 0; i < ncon; i++)
        yx[i] /= objfun_scale;
    }
    for (Int i = 0; i < nvar; i++) {
      ux[i] *= Lambda[i];
      Real li = blx[i], ui = bux[i], pxi = px[i];
      if (li > -dciInf || ui < dciInf)
        ux[i] += mu * pxi;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      ux[j] = Max(mu, minBk) * px[j];
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
    pname[8] = 0;
    problemName.assign(pname);
  }

  void Interface::analyze_J () {
    if (UseMUMPS)
      return;
    if (LJ == 0) 
      LJ = new Factor (*env);
    LJ->analyze (*J);
  }

  void Interface::cholesky_J () {
    if (UseMUMPS)
      return;
    if (LJ == 0)
      std::cerr << "analyze should be called first" << std::endl;
    LJ->factorize (*J, cholCorrection);
    cholFacs++;
    if (!env->IsPosDef()) {
      cholCorrection = 1e-12;
      LJ->factorize (*J, cholCorrection);
      cholFacs++;
    }
  }

  Real Interface::getTime () {
    timespec ts;

    clock_gettime(CLOCK_REALTIME, &ts);
    return static_cast<Real>(ts.tv_nsec)/static_cast<double>(1e9) + 
           static_cast<Real>(ts.tv_sec);
  }

#ifndef NDEBUG
  void Interface::checkInfactibility () {
    if (std::isnan(*f)) {
      throw ("f is nan");
    }
    if (std::isnan(normc)) {
      throw ("|c| is nan");
    }
    if (std::isnan(normgp)) {
      throw ("|gp| is nan");
    }
    for (Int i = 0; i < nvar; i++) {
      if ( (xx[i] > 1e10) || (xx[i] < -1e10) || (xcx[i] > 1e10) || (xcx[i] < -1e10) ) {
        Unlimited = dciTrue;
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
    if (Unlimited)
      return;
    for (Int i = 0; i < nconI; i++) {
      if ( (sx[ineqIdx[i]] > 1e10) || (sx[ineqIdx[i]] < -1e10) || (scx[ineqIdx[i]] > 1e10) || (scx[ineqIdx[i]] < -1e10) ) {
        Unlimited = dciTrue;
        break;
      }
      assert (sx[i] < cux[ineqIdx[i]]);
      assert (sx[i] > clx[ineqIdx[i]]);
      assert (scx[i] < cux[ineqIdx[i]]);
      assert (scx[i] > clx[ineqIdx[i]]);
    }
  }
#endif
  


}
