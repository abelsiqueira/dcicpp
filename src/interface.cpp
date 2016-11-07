#include "interface.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iomanip>
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
    delpointer (l_bnd);
    delpointer (u_bnd);
    delpointer (y);
    delpointer (yineq);
    delpointer (cl);
    delpointer (cu);
    delarray   (equatn);
    delarray   (ineq_index);
    delarray   (linear);
    delpointer (xc);
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
    if ( (x == 0) || (l_bnd == 0) || (u_bnd == 0) )
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
      Real yi = yx[ineq_index[i]], li = l_bndx[nvar+i], ui = u_bndx[nvar+i];
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
            << "Number of Iterations = " << iter << std::endl
            << "Elapsed Time = " << (current_time > 0 ? current_time : 0) << " s" << std::endl;
        out << "Normal iterations = " << total_normal_iteration << std::endl
            << "Average number of normal iterations = " << total_normal_iteration*1.0/iter << std::endl
            << "Maximum number of normal iterations during one iteration = "
            << max_normal_iteration << std::endl
            << "Iterations with = 1 normal iterations = " << iter_w_1_nit << std::endl
            << "Iterations with <= 1 normal iterations = " << iter_wle_1_nit << std::endl
            << "Total number of inner normal steps = " << tRest << std::endl
            << "Average number of inner normal steps = " <<
              (total_normal_iteration > 0 ? tRest/total_normal_iteration : -1) << std::endl;
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
            << "  tRest  = " << tRest << std::endl
            << "  itssmll= " << itssmll << std::endl
            << "  gap  = " << gap << std::endl
            << "----------------------" << std::endl;
        if (ncon > 0) {
          std::cout << "lambda = " << std::endl;
          y->print_more();
        }
      }
      if (display_level > 1) {
        if (nvar <= nvarshowmax) {
          out << "x: " << std::endl;
          x->print_more ();
        } else
          out << "Number of Variables is too big to show: Increase nvarshowmax" << std::endl;
        if ( (ncon > 0) && (ncon <= nconshowmax) ) {
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
    std::string flagname;
    if (filename == 0) {
      switch (exit_flag) {
        case -1:
          flagname = "assert";
          break;
        case 0:
          flagname = "convergence";
          break;
        case 1:
          flagname = "rhomax";
          break;
        case 2:
          flagname = "maxiter";
          break;
        case 3:
          flagname = "maxrest";
          break;
        case 4:
          flagname = "restfali";
          break;
        case 5:
          flagname = "shortstep";
          break;
        case 6:
          flagname = "unlimited";
          break;
        case 7:
          flagname = "timelimit";
          break;
        case 8:
          flagname = "infeasible";
          break;
        case 9:
          flagname = "nan";
          break;
        default:
          std::stringstream aux;
          aux << "Exitflag value " << exit_flag;
          throw(aux.str());
          break;
      }
      latex_name += flagname;
      if (cholesky_failed)
        latex_name += "_cholfail";
    } else
      latex_name = filename;
    file.open (latex_name.c_str(), std::ios_base::app);

    file << problemName << " & "
         << nvar << " & "
         << ncon << " & ";
    if (table_print_level > 0) {
      file << std::scientific << std::setprecision(8) << *f << " & "
           << std::scientific << std::setprecision(8) << normgp << " & "
           << std::scientific << std::setprecision(8) << normc << " & ";
    }
    file << iter << " & "
         << (current_time > 0 ? current_time : 0) << " & "
         << ((ncon > 0) ? "con" : "unc")
         << "\\\\ \\hline\n";

    file.close ();

    if (table_print_level > 1) {
      size_t endpos = problemName.find_last_not_of(" \t");
      file.open((problemName.substr(0,endpos+1)+".tab").c_str(),
          std::ios_base::out | std::ios_base::trunc);
      file << problemName << " "
        << flagname << " "
        << std::scientific << std::setprecision(6)
          << (current_time > 0 ? current_time : 1e-6) << " "
        << std::scientific << std::setprecision(6) << *f << " "
        << std::scientific << std::setprecision(6) << normgp << " "
        << std::scientific << std::setprecision(6) << normc << " "
        << (cholesky_failed ? "cholfail" : "cholok")
        << std::endl;

      file.close ();
    }
  }

  void Interface::unc_setup (Int n, Real * x, Real * bl, Real * bu) {
    set_nvar(n);
    nconE = nconI = 0;
    this->x = new Vector (*env, n+nconI);
    xx = this->x->get_doublex();
    for (Int i = 0; i < n; i++)
      xx[i] = x[i];
    xc = new Vector (*(this->x));
    xcx = xc->get_doublex();
    l_bnd = new Vector (*env, n+nconI);
    u_bnd = new Vector (*env, n+nconI);
    l_bndx = l_bnd->get_doublex();
    u_bndx = u_bnd->get_doublex();
    for (Int i = 0; i < n; i++) {
      l_bndx[i] = bl[i];
      u_bndx[i] = bu[i];
    }
    for (Int i = 0; i < n; i++) {
      if (l_bndx[i] > -dciInf || u_bndx[i] < dciInf) {
        is_bounded = dciTrue;
        break;
      }
    }

    variable_scaling = new Real[n];
    for (Int i = 0; i < (Int) n; i++)
      variable_scaling[i] = 1.0;
  }

  void Interface::con_setup (Int n, Real * x, Real * bl, Real * bu,
      Int m, Real * y, Real * cl, Real * cu, Bool * equatn) {
    set_nvar(n);
    set_ncon(m);
    this->equatn = new Bool[m];
    this->y = new Vector (*env, m, y);
    yx = this->y->get_doublex();
    nconE = nconI = 0;
    for (Int i = 0; i < m; i++) {
      Bool tmp = equatn[i];
      this->equatn[i] = tmp;
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
      for (Int i = 0; i < m; i++) {
        if (this->equatn[i] == dciFalse) {
          ineq_index[numI] = i;
          numI++;
        }
      }
    }
    this->x = new Vector (*env, n+nconI);
    xx = this->x->get_doublex();
    for (Int i = 0; i < n; i++)
      xx[i] = x[i];
    xc = new Vector (*(this->x));
    xcx = xc->get_doublex();
    l_bnd = new Vector (*env, n+nconI);
    u_bnd = new Vector (*env, n+nconI);
    l_bndx = l_bnd->get_doublex();
    u_bndx = u_bnd->get_doublex();
    for (Int i = 0; i < n; i++) {
      l_bndx[i] = bl[i];
      u_bndx[i] = bu[i];
    }
    for (Int i = 0; i < n; i++) {
      if (l_bndx[i] > -dciInf || u_bndx[i] < dciInf) {
        is_bounded = dciTrue;
        break;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      l_bndx[nvar+i] = cl[ineq_index[i]];
      u_bndx[nvar+i] = cu[ineq_index[i]];
    }

    variable_scaling = new Real[n];
    for (Int i = 0; i < (Int) n; i++)
      variable_scaling[i] = 1.0;

    this->cl = new Vector(*env, m, cl);
    clx = this->cl->get_doublex();
    this->cu = new Vector(*env, m, cu);
    cux = this->cu->get_doublex();
  }

  void Interface::set_linear (Int n, Bool * V) {
    linear = new Bool[n];
    set_ncon(n);
    nconL = nconNL = 0;
    for (Int i = 0; i < n; i++) {
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
      (*ufn) (&this->cuter_status, &nvar, xx, f);
    } else {
      (*cfn) (&this->cuter_status, &nvar, &ncon, xx, f, cx);
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
      *f -= mu*calcPen ();
      *f /= objective_scaling;
      Int numI = nvar;
      for (Int i = 0; i < ncon; i++) {
        if (equatn[i] == dciFalse) {
          cx[i] -= xx[numI];
          numI++;
        }
      }
    }
  }

  void Interface::call_fn_xc () {
    if (ncon == 0) {
      (*ufn) (&this->cuter_status, &nvar, xcx, fxc);
    } else {
      (*cfn) (&this->cuter_status, &nvar, &ncon, xcx, fxc, cx);
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
      *fxc -= mu*calcPen_xc ();
      *fxc /= objective_scaling;
      Int numI = nvar;
      for (Int i = 0; i < ncon; i++) {
        if (equatn[i] == dciFalse) {
          cx[i] -= xcx[numI];
          numI++;
        }
      }
    }
  }

  void Interface::call_ofg (Bool grad) {
    call_ofg(xx, grad);
  }

  void Interface::call_ofg (pReal px, Bool grad) {
    if (ncon == 0)
      (*uofg) (&this->cuter_status, &nvar, px, f, gx, &grad);
    else
      (*cofg) (&this->cuter_status, &nvar, px, f, gx, &grad);
    if (*f > dciInf)
      *f = dciInf;
    else if (*f < -dciInf)
      *f = -dciInf;
    if (running) {
      Real val = 0.0;
      if (grad) {
        for (Int i = nvar; i < nvar + nconI; i++)
          gx[i] = 0;
      }
      for (Int i = 0; i < nvar+nconI; i++) {
        Real xi = px[i], li = l_bndx[i], ui = u_bndx[i];
        if (li - ui > - dciTiny)
          gx[i] = 0;
        else if (li > -dciInf || ui < dciInf) {
          if ( (xi - li) < (ui - xi) ) {
            if (grad) {
              gx[i] *= (xi - li);
              gx[i] -= mu/objective_scaling;
            }
            val += log (xi - li);
          } else {
            if (grad) {
              gx[i] *= (ui - xi);
              gx[i] += mu/objective_scaling;
            }
            val += log (ui - xi);
          }
        }
        if (grad && i < nvar)
          gx[i] *= variable_scaling[i];
      }
      *f -= mu*val;
      if (objective_scaling != 1) {
        *f /= objective_scaling;
        if (grad) {
          for (Int i = 0; i < nvar; i++)
            gx[i] /= objective_scaling;
        }
      }
    }
  }

  void Interface::call_ofg_xc (Bool grad) {
    Real tmp = *f;
    call_ofg (xcx, grad);
    *fxc = *f;
    *f = tmp;
  }

  void Interface::call_ccfsg (Bool grad, Bool scale) {
    updateScaling_x();
    static bool createdvariable_scaling = false;
    pInt nnzj = new Int(0);
    (*ccfsg) (&this->cuter_status, &nvar, &ncon, xx, cx, nnzj, &amax, Jx, Jvar, Jfun, &grad);
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
        Int j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = xx[j];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -scaling_matrix[j];
            else
              Jx[*nnzj] = -1;
            Jvar[*nnzj] = j;
            Jfun[*nnzj] = i;
            (*nnzj)++;
            j++;
          }
        }
      }

      if (!createdvariable_scaling) {
        createdvariable_scaling = true;
        if (use_variable_scaling) {
          for (Int k = 0; k < *nnzj; k++){
            Int i = Jvar[k];
            if (i < nvar)
              variable_scaling[i] = Min(Max(1.0/variable_scaling[i], Jx[k]),
                  max_variable_scaling);
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
        Int numI = nvar;
        for (Int i = 0; i < ncon; i++) {
          if (equatn[i] == dciFalse) {
            Real si = xx[numI];
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
    (*ccfsg) (&this->cuter_status, &nvar, &ncon, xcx, cx, nnzj, &amax, Jx, Jvar, Jfun, &grad);
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
        Int j = nvar;
        for (Int i = 0; i < ncon; i++) { //Lines
          if (equatn[i] == dciFalse) {
            Real si = xcx[j];
            cx[i] -= si;
            if (scale)
              Jx[*nnzj] = -scaling_matrix[j];
            else
              Jx[*nnzj] = -1;
            Jvar[*nnzj] = j;
            Jfun[*nnzj] = i;
            (*nnzj)++;
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
        Int numI = nvar;
        for (Int i = 0; i < ncon; i++) {
          if (equatn[i] == dciFalse) {
            Real si = xcx[numI];
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
      if (l_bndx[i] > u_bndx[i] - dciTiny)
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
      (*uprod) (&this->cuter_status, &nvar, &gotder, xx, ppx, ux);
    else
      (*cprod) (&this->cuter_status, &nvar, &ncon, &gotder, xx, yx, ppx, ux);
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
      if (l_bndx[i] > u_bndx[i] - dciTiny) {
        ux[i] = 0;
      } else {
        ux[i] *= scaling_matrix[i];
        Real li = l_bndx[i], ui = u_bndx[i], pxi = px[i];
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
      if (l_bndx[i] > u_bndx[i] - dciTiny)
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
      (*uprod) (&this->cuter_status, &nvar, &gotder, xcx, ppx, ux);
    else
      (*cprod) (&this->cuter_status, &nvar, &ncon, &gotder, xcx, yx, ppx, ux);
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
      if (l_bndx[i] > u_bndx[i] - dciTiny) {
        ux[i] = 0;
      } else {
          ux[i] *= scaling_matrix[i];
        Real li = l_bndx[i], ui = u_bndx[i], pxi = px[i];
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
      (*unames) (&this->cuter_status, &nvar, pname, vnames);
    else {
      char gnames[10*ncon];
      (*cnames) (&this->cuter_status, &nvar, &ncon, pname, vnames, gnames);
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
      if (cholesky_correction < cholesky_base_correction)
        cholesky_correction = cholesky_base_correction;
      else
        cholesky_correction *= chol_correction_increase;
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
    for (Int i = 0; i < nvar+nconI; i++) {
      if ( (xx[i] > 1e10) || (xx[i] < -1e10) || (xcx[i] > 1e10) || (xcx[i] < -1e10) ) {
        is_unlimited = dciTrue;
        break;
      }
      if (l_bndx[i] - u_bndx[i] > - dciTiny)
        continue;

      if (xx[i] >= u_bndx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xx[i] < u_bndx[i]);

      if (xx[i] <= l_bndx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xx[i] > l_bndx[i]);

      if (xcx[i] >= u_bndx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xcx[i] < u_bndx[i]);

      if (xcx[i] <= l_bndx[i]) {
        std::cout << "ERROR" << std::endl;
        GDBSTOP ();
      }
      assert (xcx[i] > l_bndx[i]);
    }
    if (is_unlimited)
      return;
  }
#endif

}
