#ifndef dci_interface_h
#define dci_interface_h

#include "triplet.h"
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

namespace DCI {
  using namespace base_matrices;
  class Interface {
    public:
      // Constructors and Destructor
      Interface ();
      ~Interface ();

      // Problem functions
      int start ();
      int solve ();
      void show (std::ostream & = std::cout);
      void printLatex (char * = 0) const;
      void Cuter () { StartAtOne = dciTrue; };

      // Sine quibus non
      void set_x (size_t n, Real * V);
      void set_x (Vector &);
      void set_sol (size_t n, Real * V);
      void set_sol (Vector &);
      void set_lambda (size_t n, Real * V);
      void set_lambda (Vector &);
      void set_bl (size_t n, Real * V);
      void set_bl (Vector &);
      void set_bu (size_t n, Real * V);
      void set_bu (Vector &);
      void set_cl (size_t n, Real * V);
      void set_cl (Vector &);
      void set_cu (size_t n, Real * V);
      void set_cu (Vector &);
      void set_equatn (size_t n, Bool *V);
      void set_linear (size_t n, Bool *V);

      Real * get_x () const { return xx; }
      Real   get_f () const { return *f; }

      // Internal set nvar, ncon, nmax, mmax
      void set_nvar (Int n) { nvar = n; };
      void set_ncon (Int m) { ncon = m; };
      void set_nmax (Int n) { nmax = n; };
      void set_mmax (Int m) { mmax = m; };
      void set_amax (Int n) { amax = n; };

      // Function Assignment
      // Unconstrained
      void set_ufn (pufn p) { ufn = p; }
      void set_uofg (puofg p) { uofg = p; };
      void set_uprod (puprod p) { uprod = p; };
      void set_unames (punames p) { unames = p; };
      // Constrained
      void set_cfn (pcfn p) { cfn = p; };
      void set_cofg (pcofg p) { cofg = p; };
      void set_ccfsg (pccfsg p) { ccfsg = p; };
      void set_cprod (pcprod p) { cprod = p; };
      void set_ccifg (pccifg p) { ccifg = p; };
      void set_cnames (pcnames p) { cnames = p; };

    protected:
      //MUMPS Section
      DMUMPS_STRUC_C id;
      Real *mumps_matrix, *mumps_rhs;
      int *mumps_irn, *mumps_jcn, myid;
      void GDBSTOP ();
      // These are the pointer to the functions
      // Unconstrained
      pufn     ufn;     // f
      puofg    uofg;    // f and g
      puprod   uprod;   // H*p
      punames  unames;  // UNAMES
      // Constrained
      pcofg    cofg;    // f and g
      pcfn     cfn;     // f and c
      pccfsg   ccfsg;   // c and J
      pccifg   ccifg;   // c(i) and J(i,*)'
      pcprod   cprod;   // H(L)*p
      pcnames  cnames;  // CNAMES

      // Internal function calling
      // Used to call the functions above with the minimum needed parameters
      // since most are members.
      void call_fn ();
      void call_fn_xc ();
      void call_ofg (Bool = dciTrue);
      void call_ofg_xc (Bool = dciTrue);
      void call_prod (Bool, pReal, pReal);
      void call_prod_xc (Bool, pReal, pReal);
      void call_ccfsg (Bool = dciTrue, Bool = dciTrue);
      void call_ccfsg_xc (Bool = dciTrue, Bool = dciTrue);
      void call_names ();

      // Internal problem solving functions
      void analyze_J ();
      void cholesky_J ();
      Int vertstep ();
      Int vertSafeguard ();
      void horzstep (Real &);
      Int dcisteih (Vector &, Real &, Real &);
      Int dcitrust (Real);
      Int Porcelli (Real &);
      Int LeastSquareTrustRegion (Vector &, pReal, pReal, pReal);
      Int InteriorPointRestoration ();
      Real InteriorPointObjFun (Real, Real, Int, Int, pReal, pReal, pReal, pReal,
          pReal, pReal, pReal, pReal, pInt, pInt);
      void NAproj (Vector &, Vector &, Vector &);
      void NAprojApprox (Vector &, Vector &, Vector &);
      Int NAstep (Vector &, Vector &);
      Int dcibfgs (const Vector &, Int &);
      Int dcibfgs2 (const Vector &, Int &);
      void Hiprod (Int, Real, Real, pReal, pReal, pReal, pReal, pReal);
      Int linesearch (const Vector &, const Vector &, const Vector &, Real &, Real &, Vector &);
      Int zoom (const Vector &, const Vector &, const Vector &, Vector &, Real, Real, Real, Real &, Real &, Real, Real, Real, Real);
      Real interpolate (Real, Real, Real, Real, Real, Real, Bool);
      void DefineParameters ();
      void ReadParameters ();
      void InitialValues ();
      void Initialization ();
      Real getTime ();
#ifndef NDEBUG
      void checkInfactibility ();
#endif
      void UpdateScaling_x ();
      void UpdateScaling_xc ();
      void scale_x (Vector &);
      void scale_xc (Vector &);
      void project_bounds_x (Vector &);
      void project_bounds_xc (Vector &);
      void updyineq ();
      void updyineq_xc ();
      void update_lambda ();
      void update_mu ();
      Real calc_gap ();
      Real calc_ydif ();
      Real calc_pen ();
      Real calc_pen_xc ();
      Real penvthv (const Vector &) const;
      Bool calc_feasibilityOpt ();
      void LstSqrCG (Bool, const Vector &, Vector &, Vector &);
      void LinSysCG (Bool, const Vector &, Vector &);
      void assert (Bool);

      // These are the values at the current point x;
      Real * f, * fxc;
      Real mu, mugap, murho;
      Vector * g;
      Sparse * H;
      Triplet * Htrip;
      Vector * c;
      Sparse * J;
      Triplet * Jtrip;
      Factor * LJ;
      Vector * gp;
      Vector * yineq;

      Vector * x, * bl, * bu;
      Vector * solx;
      Vector * y, * cl, * cu;
      Vector * s; //Slacks s = -c_I(x)
      Vector * sols;
      Bool * equatn, * linear;
      Vector * xc, *sc;
      Vector * feasOpt;

      pReal xx, blx, bux;
      pReal yx, clx, cux;
      pReal sx;
      pReal gx, cx;
      pReal Jx;
      pInt Ji, Jj;
      pReal xcx, scx;
      pReal yineqx;

      Int * ineqIdx;
      Real * Lambda; // Scaling matrix
      Real * VariableScaling;

      Real Lref, Lprev, Lc, Ln, Lnew;
      Real DLV, DLH;
      Real gap, lagrgap, infacgap;

      Real normgp;
      Real normg, normc, normy, normx;
      Real normck;
      Real ngp; // normgp / (normg + 1)
      Real engp;
      
      // The environment variables
      Environment * env;
      Int nmax, mmax, amax;
      Int nvar, ncon, nconE, nconI, nconL, nconNL;
      std::string problemName;
      Bool StartAtOne;
      Bool Initialized;
      Bool Running;
      Bool Solved;
      Bool Ineq, Linear, Bounded;
      Bool Unlimited;
      Int DisplayLevel, VerboseLevel, TableLevel;
      Int ExitFlag;
      Real MaxTime, CurrentTime, StartTime;

      // Parameters and Variables
      Real c1, c2;
      Real rho, rhomax, rhomin;
      Real DeltaH, DeltaV;
      Real DeltaMax, DeltaMin, Delta0, DeltaTiny;
      Real alphaR, alphaI, alphaS;
      Real eps1, eps2, eps3;
      Real epsmu, epsgap;
      Real eta1, eta2, eta3;
      Real phi1, phi2;
      Real csic, csig;
      Real zeta1, zeta2, zeta3;
      Real kappa1, kappa2, kappa3, kappa4;
      Real thetaR, LbdMax;
      Bool GotH, first;
      Int nSteih, nRej, nSoc, nfailv, nHprod, nRest, nbfgs;
      Int tSoc, tSteih, tRej, tbfgs, tRest;
      Int HorzFlag, StepFlag, VertFlag;
      Int iter, maxit, maxitSteih, minitSteih, relitSteih;
      Int minstep, itssmll, maxrest, maxssmll;
      Int bfgsupd;
      Int cholFacs;
      Bool Aavail, gavail, LimLbd, FreshA;
      Real minBk;
      Bool UseCG, UsePorcelli;
      Bool UseMUMPS, ScaleVertical;
      Bool UseVertInteriorPoint, UseVertSafeguard, RebootOnVertFail;
      Bool PartialPenal, project_dcp, project_dn, project_bfgs;
      Bool PorcelliPenal;
      Bool trustWorstdn, trustConvexBox, penal_trust, penal_bfgs;
      Real cholCorrection;
      Bool cholFailed;
      Real MaxDiag, MinDiag;
      Real infeasible_gradient;
      Real choleskyCorrection;
      //Interior Point Variables

      Real objfun_scale, max_objfun_scale;
      Bool UseObjfunScale, UseVariableScaling;
      Int objfun_count;
  };
}

#endif
