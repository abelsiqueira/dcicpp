#include <interface.h>
//#include <cassert>
#include <algorithm>
#include <fstream>
#include <cmath>

namespace DCI {
  int Interface::solve () {
    if (!Initialized)
      return -1;

    //Linear problem. Enforce non-linear because the matrix is scaled.
    //Hence, we need to refactorize.
    Linear = dciFalse;
    Running = dciTrue;

#ifdef LOCALTEST
    Real ngpzk = 1.0;
    Real ngpzck = 1.0;
    Real eck = dciInf;
    Real eckp = dciInf;
    Real eckpp = dciInf;
    Real ek = dciInf;
    Real ekp = dciInf;
    Real ekpp = dciInf;
#endif

    Real norms = 0;
    Real cnormi;
    if (ncon > 0)
      cnormi = c->norm (0);
    else
      cnormi = 0;

#ifdef PLOT_MATLAB
      std::string plot_filename("plot_matlab");
      plot_filename += problemName;
      plot_filename += ".m";
      plot_filename.erase(remove_if(plot_filename.begin(),
            plot_filename.end(), isspace), plot_filename.end());
      std::ofstream plot_file(plot_filename.c_str());
      plot_file << "normc = " << (ncon > 0 ? c->norm() : 0) << ';' << std::endl
                << "ngp = " << ngp << ';' << std::endl;
#endif

#ifdef VERBOSE
    if (VerboseLevel > 0) {
      std::cout
          << "----------------------" << std::endl
          << "  DCI-C++ First Iter" << std::endl
          << "  Iter   = " << iter << std::endl
          << "  f(x)   = " << *f << std::endl
          << "  |c(x)| = " << c->norm () << std::endl
          << "  |gp|   = " << gp->norm () << std::endl
          << "  rho    = " << rho << std::endl
          << "----------------------" << std::endl;
    }
    if ( (nvar + ncon < 10) && (VerboseLevel > 1) ) {
      std::cout << "x = " << std::endl;
      x->print_more();
      std::cout << "xc = " << std::endl;
      xc->print_more();
      std::cout << "lambda = " << std::endl;
      y->print_more();
      if (Ineq) {
        std::cout << "s = " << std::endl;
        s->print_more();
        std::cout << "sc = " << std::endl;
        sc->print_more();
      }
    }
    GDBSTOP();
#endif

    CurrentTime = getTime() - StartTime;

    lagrgap = normgp;
    infacgap = normc;
    gap = 1;
    Real ydif = 0;

    while ( ( (cnormi > csic) || 
              ( (normgp > csig) && (ngp > csig*1e-2) ) || 
              (mu > epsmu) || 
              (gap > epsgap) || 
              (ydif > 1e-6) ) && 
            (iter <= maxit) && 
            (tRest <= maxrest) && 
            (itssmll <= maxssmll) && 
            (VertFlag == 0) && 
            (rhomax >= rhomin) && 
            (!Unlimited) && 
            (CurrentTime < MaxTime) ) {

      calc_ydif ();
      iter++;

      //Check for infactibility
#ifndef NDEBUG
      checkInfactibility();
#endif

      if (DeltaH < DeltaMin) DeltaH = DeltaMin;
      if (DeltaV < DeltaMin) DeltaV = DeltaMin;

      vertstep (); //Recalculates f, g and c

#ifdef LOCALTEST
      ngpzck = gp->norm();
      if (solx != 0) {
        eck = eckp;
        eckp = eckpp;
        eckpp = norm(*xc - *solx) + norm(*sc - *sols);
      }
#endif

      tRest += nRest;
      tbfgs += nbfgs;
//      assert ( (normc <= rho) || (VertFlag != 0) || (tRest > maxrest) );

      call_ofg_xc (dciFalse);

      if (ncon > 0)
        cnormi = c->norm (0);
      else
        cnormi = 0;

      if (ncon > 0)
        Lc = *fxc + y->dot (*c);
      else
        Lc = *fxc;
      DLV = Lc - Ln;

      if (DLV >= 0.5*(Lref - Ln))
        rhomax = rhomax/2;
      if (DLV > -0.5*DLH)
        Lref = Lc;

#ifdef PLOT_MATLAB
      plot_file << "normc(end+1) = " << (ncon > 0 ? c->norm() : 0) << ';' << std::endl
                << "ngp(end+1) = " << ngp << ';' << std::endl;
#endif

#ifdef VERBOSE
    if (VerboseLevel > 0) {
      std::cout << std::endl
          << "----------------------" << std::endl
          << "Vertical Step" << std::endl
          << "  Iter   = " << iter << std::endl
          << "  f(x)   = " << *fxc << std::endl
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
          << "  stepsize = " << (*x - *xc).norm() << std::endl
          << "----------------------" << std::endl;
    }
    if ( (nvar + ncon < 10) && (VerboseLevel > 1) ) {
      std::cout << "x = " << std::endl;
      x->print_more();
      std::cout << "xc = " << std::endl;
      xc->print_more();
      std::cout << "lambda = " << std::endl;
      y->print_more();
      if (Ineq) {
        std::cout << "s = " << std::endl;
        s->print_more();
        std::cout << "sc = " << std::endl;
        sc->print_more();
      }
    }
    GDBSTOP();
#endif
#ifdef LOCALTEST
    if ( (solx != 0) && (iter > 3) ) {
      std::cout << "eckpp/eck = " << eckpp/eck << std::endl;
    }
#endif

      if (VertFlag != 0)
        std::cout << "VertFlag = " << VertFlag << std::endl;

      CurrentTime = getTime() - StartTime;

      if ( ( (cnormi > csic) || 
             ( (normgp > csig) && 
               (ngp > csig*1e-2) ) ) && 
           (VertFlag == 0) && 
           (rhomax >= rhomin) && 
           (!Unlimited) && 
           (CurrentTime < MaxTime) ) {

        horzstep (norms); 
#ifndef NDEBUG
        checkInfactibility();
#endif

        normck = normc;

        tSoc += nSoc;
        tSteih += nSteih;
        tRej += nRej;


        if (nSteih > 0)
          engp = fabs (DLH)/( fabs((*f) - (*fxc)) + norms);
        else
          engp = ngp;

        normx = x->norm();
        if (norms < normx * minstep)
          itssmll++;
        else
          itssmll = 0;

      } else {
        *x = *xc;
        if (Ineq)
          *s = *sc;
        *f = *fxc;
      } 
      updyineq ();
      gap = calc_gap ();

#ifdef LOCALTEST
      Vector ytmp (*y);
      call_ofg(dciTrue);
      update_lambda ();
      ngpzk = gp->norm();
      *y = ytmp;
      if (solx != 0) {
        ek = ekp;
        ekp = ekpp;
        ekpp = norm(*x - *solx) + norm(*s - *sols);
      }
#endif

#ifdef PLOT_MATLAB
      plot_file << "normc(end+1) = " << (ncon > 0 ? c->norm() : 0) << ';' << std::endl
                << "ngp(end+1) = " << ngp << ';' << std::endl;
#endif
#ifdef VERBOSE
      if (VerboseLevel > 0) {
        std::cout << std::endl
            << "----------------------" << std::endl
            << "Horizontal Step" << std::endl
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
            << "  stepsize = " << (*x - *xc).norm() << std::endl
            << "----------------------" << std::endl;
      }
      if ( (nvar + ncon < 10) && (VerboseLevel > 1) ) {
        std::cout << "x = " << std::endl;
        x->print_more();
        std::cout << "xc = " << std::endl;
        xc->print_more();
        std::cout << "lambda = " << std::endl;
        y->print_more();
        if (Ineq) {
          std::cout << "s = " << std::endl;
          s->print_more();
          std::cout << "sc = " << std::endl;
          sc->print_more();
        }
      }
      GDBSTOP();
#endif
#ifdef LOCALTEST
      std::cout << "|gp(zk)|/|gp(zck)| = " << ngpzk/ngpzck << std::endl;
      if ( (solx != 0) && (iter > 3) ) {
        std::cout << "ekpp/ek = " << ekpp/ek << std::endl;
      }
#endif
      
    } //End

    if (VertFlag == 2)
      ExitFlag = 8;
    else if (VertFlag > 0)
      ExitFlag = 4;
    else if (rhomax < rhomin)
      ExitFlag = 1;
    else if (iter > maxit)
      ExitFlag = 2;
    else if (tRest > maxrest)
      ExitFlag = 3;
    else if (itssmll > maxssmll)
      ExitFlag = 5;
    else if (Unlimited)
      ExitFlag = 6;
    else if (CurrentTime >= MaxTime)
      ExitFlag = 7;
    else
      ExitFlag = 0;

    Running = dciFalse;
    Solved = dciTrue;

    //Calculating the real function value
    call_fn();
    CurrentTime = getTime() - StartTime;

    //MUMPS
    if (UseMUMPS) {
      id.job=JOB_END;
      dmumps_c(&id);
      MPI_Finalize();
    }

#ifdef PLOT_MATLAB
    plot_file.close();
#endif

    return 0;
  }
}
