#include "interface.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <fstream>
//#include <cassert>

namespace DCI {
  void error (int status, const char *file, int line, const char *message) {
#ifdef VERBOSE
     std::cerr << message << std::endl
               << "Status = " << status << std::endl
               << "File   = " << file << std::endl
               << "Line   = " << line << std::endl;
#else
     UNUSED(status);
     UNUSED(file);
     UNUSED(line);
     UNUSED(message);
#endif
}

  void Interface::assert (Bool v) {
    if (v)
      return;
    exit_flag = -1;
    error (-1, "assert", -1, "assert failed");
    GDBSTOP ();
    throw "assert failed";
  }

  void Interface::GDBSTOP () {
    return;
  }

  Interface::Interface () {
    //Parameters
    initialization ();

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
    delarray   (ineq_index);
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
    delpointer (scaling_matrix);
    delpointer (variable_scaling);
  }

  int Interface::start () {
    //Tests for null
    if ( (x == 0) || (bl == 0) || (bu == 0) )
      return -1;

    call_names ();
    start_time = getTime();

    if ( (ncon > 0) && (equatn == 0) ) {
      equatn = new Bool[ncon];
      for (Int i = 0; i < ncon; i++)
        equatn[i] = dciTrue;
      nconE = ncon;
      nconI = 0;
      has_ineq = dciFalse;
    }

    if (linear == 0) {
      linear = new Bool[ncon];
      for (Int i = 0; i < ncon; i++)
        linear[i] = dciFalse;
      nconNL = ncon;
      nconL = 0;
      is_linear = dciFalse;
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
      if (has_ineq) {
        s = new Vector (*env, nconI);
        sx = s->get_doublex();
        sc = new Vector (*env, nconI);
        scx = sc->get_doublex();
      }
      Jtrip = new Triplet (*env, mmax, nmax + nconI, amax, 0); //0 is for unsymmetric
      Jx = Jtrip->get_doublex();
      Jfun = static_cast < Int * > (Jtrip->triplet->i);
      Jvar = static_cast < Int * > (Jtrip->triplet->j);
      for (Int i = 0; i < amax; i++) {
        Jx[i] = 0;
        Jfun[i] = 0;
        Jvar[i] = 0;
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
        ps[i] = cx[ineq_index[i]];
      }
    }
#endif

    defineParameters ();
    initialValues ();
    initialized = dciTrue;

    return 0;
  }

  void Interface::show (std::ostream & out) {
    if (!initialized) {
      out << "Problem not initialized" << std::endl;
      return;
    }

    Real yoff = 0;
    for (Int i = 0; i < nconI; i++) {
      Real yi = yx[ineq_index[i]], li = clx[ineq_index[i]], ui = cux[ineq_index[i]];
      if (li > -dciInf)
        yoff += Max(yi, 0.0);
      if (ui < dciInf)
        yoff += Max(-yi, 0.0);
    }

    if (display_level > 0) {
      out << "**********************************************************" << std::endl;
      out << "Problem name: " << problemName << std::endl
          << std::endl
          << "Number of Variables: " << nvar << std::endl
          << "Number of Constraints: " << ncon << std::endl 
          << "             Equality: " << nconE << std::endl
          << "           has_inequality: " << nconI << std::endl
          << std::endl;
    }

#ifdef PRINT_MATLAB
    if (J != 0) {
      std::string matlab_filename("matrix_jacob");
      matlab_filename += problemName;
      switch (exit_flag) {
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

    if (!solved) {
      out << "Problem not solved yet" << std::endl;
    } else {
      out << "EXIT: ";
      if (exit_flag == 0)
        out << "The Algorithm has Converged" << std::endl;
      else if (exit_flag == 1)
        out << "rhomax became too short" << std::endl;
      else if (exit_flag == 2)
        out << "The maximum number of iterations was reached" << std::endl;
      else if (exit_flag == 3)
        out << "The maximum number of restorarions was reached" << std::endl;
      else if (exit_flag == 4)
        out << "The restoration has failed" << std::endl;
      else if (exit_flag == 5)
        out << "The step became too short" << std::endl;
      else if (exit_flag == 6)
        out << "The problem is unlimited" << std::endl;
      else if (exit_flag == 7)
        out << "The problem reached the time limit" << std::endl;
      else if (exit_flag == 8)
        out << "Stopped at a stationary for the infeasibility" << std::endl;

      if (cholesky_failed)
        out << "Cholesky failed" << std::endl;

      if (display_level > 0) {
        out << "f(x) = " << *f << std::endl
            << "|c(x)| = " << normc << std::endl
            << "|g(x) + J(x)'*y| = " << normgp << std::endl
            << "y offset = " << yoff << std::endl
            << "BFGS? " << ((tbfgs > 0) ? "yes" : "no") << std::endl
            << "Number of Iterations = " << iter << std::endl
            << "Elapsed Time = " << (current_time > 0 ? current_time : 0) << " s" << std::endl;
      }

      if (display_level > 2) {
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
          if (has_ineq) {
            std::cout << "s = " << std::endl;
            s->print_more();
            std::cout << "sc = " << std::endl;
            sc->print_more();
          }
        }
      }
      if (display_level > 1) {
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

  /* 
   * If (table_print_level == 0)
   * Problem name & nvar & ncon & iters & time & converged
   * If (table_print_level == 1)
   * Problem name & nvar & ncon & fvalue & normgp & normh & iters & time & converged
   */
  void Interface::printLatex (char * filename) const {
    std::ofstream file;
    std::string latex_name("latex_");
    if (filename == 0) {
      switch (exit_flag) {
        case -1:
          latex_name += "assert";
          break;
        case 0:
          latex_name += "convergence";
          break;
        case 1:
          latex_name += "rhomax";
          break;
        case 2:
          latex_name += "maxiter";
          break;
        case 3:
          latex_name += "maxrest";
          break;
        case 4:
          latex_name += "restfali";
          break;
        case 5:
          latex_name += "shortstep";
          break;
        case 6:
          latex_name += "unlimited";
          break;
        case 7:
          latex_name += "timelimit";
          break;
        case 8:
          latex_name += "infeasible";
          break;
        case 9:
          latex_name += "nan";
          break;
        default:
          std::stringstream aux;
          aux << "Exitflag value " << exit_flag;
          throw(aux.str());
          break;
      }
      if (cholesky_failed)
        latex_name += "_cholfail";
    } else
      latex_name = filename;
    file.open (latex_name.c_str(), std::ios_base::app);
    
    file << problemName << " & "
         << nvar << " & "
         << ncon << " & ";
    if (table_print_level > 0) {
      file << *f << " & "
           << normgp << " & "
           << normc << " & ";
    }
    file << iter << " & "
         << (current_time > 0 ? current_time : 0) << " & "
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
    variable_scaling = new Real[n];
    for (Int i = 0; i < (Int)n; i++)
      variable_scaling[i] = 1.0;
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
        is_bounded = dciTrue;
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
        is_bounded = dciTrue;
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
      has_ineq = dciTrue;
      is_bounded = dciTrue;
      ineq_index = new Int[nconI];
      Int numI = 0;
      for (size_t i = 0; i < n; i++) {
        if (equatn[i] == dciFalse) {
          ineq_index[numI] = i;
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
      is_linear = dciTrue;
  }

  void Interface::call_fn () {
    if (ncon == 0) {
      (*ufn) (&nvar, xx, f);
    } else {
      (*cfn) (&nvar, &ncon, xx, f, &mmax, cx);
      for (Int i = 0; i < ncon; i++) {
        if (cx[i] > dciInf)
          cx[i] = dciInf;
        else if (cx[i] < -dciInf)
          cx[i] = -dciInf;
      }
    }
    if (ncon > 0 && use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        cx[i] /= constraint_scaling[i];
    }
    if (*f > dciInf)
      *f = dciInf;
    else if (*f < -dciInf)
      *f = -dciInf;
    if (running) {
      *f /= objective_scaling;
      *f -= mu*calcPen ();
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
      (*ufn) (&nvar, xcx, fxc);
    } else {
      (*cfn) (&nvar, &ncon, xcx, fxc, &mmax, cx);
      for (Int i = 0; i < ncon; i++) {
        if (cx[i] > dciInf)
          cx[i] = dciInf;
        else if (cx[i] < -dciInf)
          cx[i] = -dciInf;
      }
    }
    if (ncon > 0 && use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        cx[i] /= constraint_scaling[i];
    }
    if (*fxc > dciInf)
      *fxc = dciInf;
    else if (*fxc < -dciInf)
      *fxc = -dciInf;
    if (running) {
      *fxc /= objective_scaling;
      *fxc -= mu*calcPen_xc ();
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
    if (*f > dciInf)
      *f = dciInf;
    else if (*f < -dciInf)
      *f = -dciInf;
    if (running) {
      Real val = 0.0;
      if (objective_scaling != 1) {
        *f /= objective_scaling;
        if (grad) {
          for (Int i = 0; i < nvar; i++)
            gx[i] /= objective_scaling;
        }
      }
      for (Int i = 0; i < nvar; i++) {
        Real xi = xx[i], bli = blx[i], bui = bux[i];
        if (bli - bui > - dciTiny) 
          gx[i] = 0;
        else if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (partial_penalization) {
            if ( (xi - bli) < (bui - xi) ) {
              if (grad) {
                gx[i] *= (xi - bli);
                gx[i] -= mu/objective_scaling;
              }
              val += log (xi - bli);
            } else {
              if (grad) {
                gx[i] *= (bui - xi);
                gx[i] += mu/objective_scaling;
              }
              val += log (bui - xi);
            }
          } else {
            if (grad) {
              gx[i] *= (xi - bli) * (bui - xi);
              gx[i] += mu * (2*xi - bli - bui)/objective_scaling;
            }
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          if (grad) {
            gx[i] *= (bui - xi);
            gx[i] += mu/objective_scaling;
          }
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          if (grad) {
            gx[i] *= (xi - bli);
            gx[i] -= mu/objective_scaling;
          }
          val += log (xi - bli);
        }
        if (grad)
          gx[i] *= variable_scaling[i];
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = sx[i], cli = clx[ineq_index[i]], cui = cux[ineq_index[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (partial_penalization) {
            if ( (si - cli) < (cui - si) ) {
              if (grad)
                gx[j] = -mu/objective_scaling;
              val += log (si - cli);
            } else {
              if (grad)
                gx[j] = mu/objective_scaling;
              val += log (cui - si);
            }
          } else {
            if (grad)
              gx[j] = mu * (2*si - cli - cui)/objective_scaling;
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          if (grad)
            gx[j] = mu/objective_scaling;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          if (grad)
            gx[j] = -mu/objective_scaling;
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
    if (*fxc > dciInf)
      *fxc = dciInf;
    else if (*fxc < -dciInf)
      *fxc = -dciInf;
    if (running) {
      Real val = 0.0;
      if (objective_scaling != 1) {
        *fxc /= objective_scaling;
        if (grad) {
          for (Int i = 0; i < nvar; i++)
            gx[i] /= objective_scaling;
        }
      }
      for (Int i = 0; i < nvar; i++) {
        Real xi = xcx[i], bli = blx[i], bui = bux[i];
        if (bli - bui > - dciTiny) 
          gx[i] = 0;
        else if ( (bli > -dciInf) && (bui < dciInf) ) {
          if (partial_penalization) {
            if ( (xi - bli) < (bui - xi) ) {
              if (grad) {
                gx[i] *= (xi - bli);
                gx[i] -= mu/objective_scaling;
              }
              val += log (xi - bli);
            } else {
              if (grad) {
                gx[i] *= (bui - xi);
                gx[i] += mu/objective_scaling;
              }
              val += log (bui - xi);
            }
          } else {
            if (grad) {
              gx[i] *= (xi - bli) * (bui - xi);
              gx[i] += mu * (2*xi - bli - bui)/objective_scaling;
            }
            val += log (xi - bli) + log (bui - xi);
          }
        } else if ( (bli <= -dciInf) && (bui < dciInf) ) {
          if (grad) {
            gx[i] *= (bui - xi);
            gx[i] += mu/objective_scaling;
          }
          val += log (bui - xi);
        } else if ( (bli > -dciInf) && (bui >= dciInf) ) {
          if (grad) {
            gx[i] *= (xi - bli);
            gx[i] -= mu/objective_scaling;
          }
          val += log (xi - bli);
        }
        if (grad)
          gx[i] *= variable_scaling[i];
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real si = scx[i], cli = clx[ineq_index[i]], cui = cux[ineq_index[i]];
        if ( (cli > -dciInf) && (cui < dciInf) ) {
          if (partial_penalization) {
            if ( (si - cli) < (cui - si) ) {
              if (grad)
                gx[j] = -mu/objective_scaling;
              val += log (si - cli);
            } else {
              if (grad)
                gx[j] = mu/objective_scaling;
              val += log (cui - si);
            }
          } else {
            if (grad)
              gx[j] = mu * (2*si - cli - cui)/objective_scaling;
            val += log (si - cli) + log (cui - si);
          }
        } else if ( (cli <= -dciInf) && (cui < dciInf) ) {
          if (grad)
            gx[j] = mu/objective_scaling;
          val += log (cui - si);
        } else if ( (cli > -dciInf) && (cui >= dciInf) ) {
          if (grad)
            gx[j] = -mu/objective_scaling;
          val += log (si - cli);
        }
      }
      *fxc -= mu*val;
    }
  }

  void Interface::call_ccfsg (Bool grad, Bool scale) {
    updateScaling_x();
    static bool createdvariable_scaling = false;
    pInt nnzj = new Int(0);
    (*ccfsg) (&nvar, &ncon, xx, &mmax, cx, nnzj, &amax, Jx, Jvar, Jfun, &grad);
    for (Int i = 0; i < ncon; i++) {
      if (cx[i] > dciInf)
        cx[i] = dciInf;
      else if (cx[i] < -dciInf)
        cx[i] = -dciInf;
    }
    for (Int i = 0; i < *nnzj; i++) {
      if (Jx[i] > dciInf)
        Jx[i] = dciInf;
      else if (Jx[i] < -dciInf)
        Jx[i] = -dciInf;
    }
    if (grad == dciTrue) {
      if (start_at_one) {
        for (Int i = 0; i < *nnzj; i++) {
          Jvar[i]--;
          Jfun[i]--;
        }
      }
      // Fixed variables fix
      for (Int k = 0; k < *nnzj; k++) {
        Real bli = blx[Jvar[k]], bui = bux[Jvar[k]];
        if (bli - bui > - dciTiny) {
          Jx[k]   = Jx[*nnzj-1];
          Jvar[k] = Jvar[*nnzj-1];
          Jfun[k] = Jfun[*nnzj-1];
          (*nnzj)--;
        }
      }

      if (running) {
        if (ncon > 0 && use_constraint_scaling) {
          for (Int i = 0; i < ncon; i++)
            cx[i] /= constraint_scaling[i];
          for (Int k = 0; k < *nnzj; k++) {
            Int i = Jfun[k];
            Jx[k] /= constraint_scaling[i];
          }
        }
        if (scale) { 
          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jvar[k];
            Jx[k] *= scaling_matrix[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = sx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -scaling_matrix[j];
            else
              Jx[*nnzj] = -1;
            Jvar[*nnzj] = j;
            Jfun[*nnzj] = i;
            (*nnzj)++;
            numI++;
            j++;
          }
        }
      }

      if (!createdvariable_scaling) {
        createdvariable_scaling = true;
        if (use_variable_scaling) {
          for (Int k = 0; k < *nnzj; k++){
            Int i = Jvar[k];
            variable_scaling[i] = Max(variable_scaling[i], Jx[k]);
          }
          updateScaling_x();
        }
        for (Int i = 0; i < nvar; i++)
          variable_scaling[i] = 1.0/variable_scaling[i];
      }

      *(Jtrip->get_pnnz()) = *nnzj;
      delpointer (J);
      J = new Sparse (*env);
      J->triplet_to_sparse (*Jtrip, amax);
    } else {
      if (ncon > 0 && use_constraint_scaling) {
        for (Int i = 0; i < ncon; i++)
          cx[i] /= constraint_scaling[i];
      }
      if (running) {
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
    updateScaling_xc();
    pInt nnzj = new Int(0);
    (*ccfsg) (&nvar, &ncon, xcx, &mmax, cx, nnzj, &amax, Jx, Jvar, Jfun, &grad);
    for (Int i = 0; i < ncon; i++) {
      if (cx[i] > dciInf)
        cx[i] = dciInf;
      else if (cx[i] < -dciInf)
        cx[i] = -dciInf;
    }
    for (Int i = 0; i < *nnzj; i++) {
      if (Jx[i] > dciInf)
        Jx[i] = dciInf;
      else if (Jx[i] < -dciInf)
        Jx[i] = -dciInf;
    }
    if (grad == dciTrue) {
      if (start_at_one) {
        for (Int i = 0; i < *nnzj; i++) {
          Jvar[i]--;
          Jfun[i]--;
        }
      }
      // Fixed variables fix
      for (Int k = 0; k < *nnzj; k++) {
        Real bli = blx[Jvar[k]], bui = bux[Jvar[k]];
        if (bli - bui > - dciTiny) {
          Jx[k]   = Jx[*nnzj-1];
          Jvar[k] = Jvar[*nnzj-1];
          Jfun[k] = Jfun[*nnzj-1];
          (*nnzj)--;
        }
      }

      if (running) {
        if (ncon > 0 && use_constraint_scaling) {
          for (Int i = 0; i < ncon; i++)
            cx[i] /= constraint_scaling[i];
          for (Int k = 0; k < *nnzj; k++) {
            Int i = Jfun[k];
            Jx[k] /= constraint_scaling[i];
          }
        }
        if (scale) {
          for (Int k = 0; k < *nnzj; k++) {
            Int j = Jvar[k];
            Jx[k] *= scaling_matrix[j];
          }
        }
        Int numI = 0, j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = scx[numI];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -scaling_matrix[j];
            else
              Jx[*nnzj] = -1;
            Jvar[*nnzj] = j;
            Jfun[*nnzj] = i;
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
      if (ncon > 0 && use_constraint_scaling) {
        for (Int i = 0; i < ncon; i++)
          cx[i] /= constraint_scaling[i];
      }
      if (running) {
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
    updateScaling_x();
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > bux[i] - dciTiny)
        ppx[i] = 0.0;
      else
        ppx[i] = scaling_matrix[i] * px[i];
    }
    if (objective_scaling != 1) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= objective_scaling;
    }
    if (use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        yx[i] /= constraint_scaling[i];
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xx, &mmax, yx, ppx, ux);
    if (objective_scaling != 1) {
      for (Int i = 0; i < nvar; i++)
        ux[i] /= objective_scaling;
      for (Int i = 0; i < ncon; i++)
        yx[i] /= objective_scaling;
    }
    if (use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= constraint_scaling[i];
    }
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > bux[i] - dciTiny) {
        ux[i] = 0;
      } else {
        ux[i] *= scaling_matrix[i];
        Real li = blx[i], ui = bux[i], pxi = px[i];
        if (li > -dciInf || ui < dciInf)
          ux[i] += mu * pxi;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      ux[j] = Max(mu, minBk) * px[j];
    }
  }

  void Interface::call_prod_xc (Bool gotder, pReal px, pReal ux) {
    updateScaling_xc();
    Real ppx[nvar];
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > bux[i] - dciTiny)
        ppx[i] = 0;
      else
        ppx[i] = scaling_matrix[i] * px[i];
    }
    if (objective_scaling != 1) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= objective_scaling;
    }
    if (use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        yx[i] /= constraint_scaling[i];
    }
    if (ncon == 0)
      (*uprod) (&nvar, &gotder, xcx, ppx, ux);
    else
      (*cprod) (&nvar, &ncon, &gotder, xcx, &mmax, yx, ppx, ux);
    if (objective_scaling != 1) {
      for (Int i = 0; i < nvar; i++)
        ux[i] /= objective_scaling;
      for (Int i = 0; i < ncon; i++)
        yx[i] /= objective_scaling;
    }
    if (use_constraint_scaling) {
      for (Int i = 0; i < ncon; i++)
        yx[i] *= constraint_scaling[i];
    }
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > bux[i] - dciTiny) {
        ux[i] = 0;
      } else {
          ux[i] *= scaling_matrix[i];
        Real li = blx[i], ui = bux[i], pxi = px[i];
        if (li > -dciInf || ui < dciInf)
          ux[i] += mu * pxi;
      }
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

  void Interface::analyzeJacobian () {
    if (LJ == 0) 
      LJ = new Factor (*env);
    LJ->analyze (*J);
  }

  void Interface::cholesky () {
    if (LJ == 0)
      std::cerr << "analyze should be called first" << std::endl;
    LJ->factorize (*J, cholesky_correction);
    cholFacs++;
    if (!env->IsPosDef()) {
      cholesky_failed = dciTrue;
      cholesky_correction = 1e-12;
      LJ->factorize (*J, cholesky_correction);
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
    try {
      if (std::isnan(*f)) {
        throw ("f is nan");
      }
      if (std::isnan(normc)) {
        throw ("|c| is nan");
      }
      if (std::isnan(normgp)) {
        throw ("|gp| is nan");
      }
      if (std::isnan(x->norm())) {
        throw ("|x| is nan");
      }
      if (std::isnan(xc->norm())) {
        throw ("|xc| is nan");
      }
    } catch (const char * ex) {
      throw(ex);
    }
    for (Int i = 0; i < nvar; i++) {
      if ( (xx[i] > 1e10) || (xx[i] < -1e10) || (xcx[i] > 1e10) || (xcx[i] < -1e10) ) {
        is_unlimited = dciTrue;
        break;
      }
      if (blx[i] - bux[i] > - dciTiny)
        continue;

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
    if (is_unlimited)
      return;
    for (Int i = 0; i < nconI; i++) {
      if ( (sx[ineq_index[i]] > 1e10) || (sx[ineq_index[i]] < -1e10) || (scx[ineq_index[i]] > 1e10) || (scx[ineq_index[i]] < -1e10) ) {
        is_unlimited = dciTrue;
        break;
      }
      assert (sx[i] < cux[ineq_index[i]]);
      assert (sx[i] > clx[ineq_index[i]]);
      assert (scx[i] < cux[ineq_index[i]]);
      assert (scx[i] > clx[ineq_index[i]]);
    }
  }
#endif

}
