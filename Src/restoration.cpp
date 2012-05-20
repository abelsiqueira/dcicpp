//#define RESTORATIONPRINT

#include "interface.h"
#include <cmath>
/* Os limites tem que ser atualizados para os pontos atuais
 *
 * l <= d + alpha*deltaD <= u */

/* Interior Point Restoration
 *
 * This routine restores the iteration to the cylinder.
 * This computes approximately a solution to
 *
 * min 0.5*||A*d + h||^2
 * s.t l <= C*d <= u
 *     ||d|| <= DeltaV
 *
 * where A is the jacobian, or scaled jacobian, and
 * C is the identity, or the Scale matrix, respectively.
 * The solution is computed by applying a newton step
 * to the KKT equations with sigma*mu:
 *
 *      | A'*(A*d + h) + C'*tu - C'*tl |
 * F =  | C*d - y                      |
 *      | (U - Y)*tu - sigma*mu*e      |
 *      | (L - Y)*tl - sigma*mu*e      |
 * */

/* Important informations:
 * A is ncon x (nvar + nconI) matrix 
 * C is diagonal. */


namespace DCI {
  Real Interface::InteriorPointObjFun (Real alphaPrimal, Real alphaDual, 
      Int numUpper, Int numLower, pReal oldrdx, pReal oldrpx, pReal dsvx, 
      pReal Cx, pReal dStep, pReal slackStep, pReal upperStep, pReal lowerStep, 
      pInt upperIndex, pInt lowerIndex) {

    Real rdx[nvar + nconI], rpx[nvar + nconI];

    for (int i = 0; i < nvar + nconI; i++) {
      rpx[i] = oldrpx[i] + alphaPrimal*(Cx[i]*dStep[i] - slackStep[i]);
      rdx[i] = oldrdx[i] + alphaPrimal*dsvx[i];
    }
    for (int i = 0; i < numUpper; i++) {
      int j = upperIndex[i];
      rdx[j] += alphaDual*Cx[j]*upperStep[i];
    }
    for (int i = 0; i < numLower; i++) {
      int j = lowerIndex[i];
      rdx[j] -= alphaDual*Cx[j]*lowerStep[i];
    }

    Real f = 0.0;

    for (int i = 0; i < nvar + nconI; i++) {
      f += pow(rpx[i], 2) + pow(rdx[i], 2);
    }

    return f;
  }

  Int Interface::InteriorPointRestoration () {
    Vector d(*env);
    Real *slack, *multUpper, *multLower;
    Real upper[nvar + nconI], lower[nvar + nconI]; //Bounds for C*d.
    Vector matrixC(*env);
    matrixC.reset(nvar + nconI, 1.0);
    Int iter = 0, maxIterIPRestoration = 10;
    Vector normGrad(*env);
    Real one[2] = {1,0}, zero[2] = {0,0};
    normGrad.sdmult(*J, 1, one, zero, *c);
    pReal normgx = 0;
    normgx = normGrad.get_doublex();
    Real oldnormc = normc;

    Vector oldxc(*xc), oldsc(*sc);
    pReal oldxcx = oldxc.get_doublex();
    pReal oldscx = oldsc.get_doublex();

    d.sdmult(*J, 0, one, zero, normGrad);
    Real alphacp = normGrad.norm()/d.norm();
    alphacp *= alphacp;
    alphacp = Min(alphacp, DeltaV/normGrad.norm());
    d.scale(normGrad, -alphacp);
//    d.reset(nvar + nconI, 0.0);
    pReal dx = d.get_doublex();

    Int numUpper = 0, numLower = 0;
    Int upperIndex[nvar + nconI], lowerIndex[nvar + nconI];

    if (ScaleVertical)
      scale_xc(matrixC);

#ifdef RESTORATIONPRINT
    std::cout << "xc = " << std::endl;
    xc->print_more();
    std::cout << "sc = " << std::endl;
    sc->print_more();
    std::cout << "|c| = " << c->norm() << std::endl;
#endif
    checkInfactibility ();

    pReal Cx = matrixC.get_doublex();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], li = blx[i], ui = bux[i];
      if (li > -dciInf) {
        lower[numLower] = Max( (li - zi) * (1 - 1e-1), -DeltaV);
        lowerIndex[numLower++] = i;
      } else {
        lower[numLower] = -DeltaV;
        lowerIndex[numLower++] = i;
      }
      if (ui < dciInf) {
        upper[numUpper] = Min( (ui - zi) * (1 - 1e-1), DeltaV);
        upperIndex[numUpper++] = i;
      } else {
        upper[numUpper] = DeltaV;
        upperIndex[numUpper++] = i;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      if (li > -dciInf) {
        lower[numLower] = (li - zi) * (1 - 1e-1);
        lowerIndex[numLower++] = j;
      }
      if (ui < dciInf) {
        upper[numUpper] = (ui - zi) * (1 - 1e-1);
        upperIndex[numUpper++] = j;
      }
    }

    Real innerMu = 1.0;
    multUpper = new Real[numUpper];
    multLower = new Real[numLower];
    slack     = new Real[nvar + nconI];
    Real rhsClone[nvar + nconI + numUpper + numLower];

    for (int i = 0; i < nvar + nconI; i++)
      slack[i] = Cx[i]*dx[i];
    for (int i = 0; i < numUpper; i++) {
      int j = upperIndex[i];
      multUpper[i] = 1.0;
      if (slack[j] >= upper[i])
        slack[j] = upper[i]/2;
    }
    for (int i = 0; i < numLower; i++) {
      int j = lowerIndex[i];
      multLower[i] = 1.0;
      if (slack[j] <= lower[i])
        slack[j] = lower[i]/2;
    }
    for (int i = 0; i < nvar + nconI; i++)
      dx[i] = slack[i]/Cx[i];

    Real normSlack = 0.0;
    for (int i = 0; i < nvar + nconI; i++)
      normSlack += pow(slack[i], 2);

    DMUMPS_STRUC_C vertID;
    int intPointMyID;

    MPI_Comm_rank(MPI_COMM_WORLD, &intPointMyID);
    vertID.job = JOB_INIT;
    vertID.par = 1;
    vertID.sym = 2;
    vertID.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&vertID);
    vertID.n = nvar + nconI + numUpper + numLower;

    Sparse G(*env); // G will be A'*A;
    G.transpose(*J); // G = A';
    Int fset[ncon];
    for (Int i = 0; i < ncon; i++)
      fset[i] = i;
    G.aat (G, fset, ncon); // G = G*G' = A'*A;
    Triplet GT(*env);
    GT.sparse_to_triplet(G);
    //Create initial matrix
    vertID.nz  = (int)(G.nnz()) + 2*(numUpper + numLower);

    vertID.irn = new int[vertID.nz];
    vertID.jcn = new int[vertID.nz];
    vertID.a   = new double[vertID.nz];

    for (int i = 0; i < numUpper; i++) {
      vertID.irn[2*i] = nvar + nconI + i + 1;
      vertID.jcn[2*i] = i + 1;
      vertID.a[2*i]   = Cx[upperIndex[i]];

      vertID.irn[2*i + 1] = nvar + ncon + i + 1;
      vertID.jcn[2*i + 1] = nvar + ncon + i + 1;
      vertID.a[2*i + 1]   = -(upper[i] - slack[upperIndex[i]])/multUpper[i];
    }
    for (int i = 0; i < numLower; i++) {
      vertID.irn[2*i + 2*numUpper] = nvar + nconI + numUpper + i + 1;
      vertID.jcn[2*i + 2*numUpper] = i + 1;
      vertID.a[2*i + 2*numUpper]   = -Cx[lowerIndex[i]];

      vertID.irn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
      vertID.jcn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
      vertID.a[2*i + 1 + 2*numUpper]   = (lower[i] - slack[lowerIndex[i]])/multLower[i];
    }
    int k = 2*numUpper + 2*numLower;
    int gtnnz = (*GT.get_pnnz());
    long int *pirn = GT.get_linti(), *pjcn = GT.get_lintj();
    double *pa = GT.get_doublex();
    bool diagonalVisited[nvar + nconI];
    for (int i = 0; i < nvar + nconI; i++)
      diagonalVisited[i] = false;
    for (int i = 0; i < gtnnz; i++) {
      if (pirn[i] < pjcn[i])
        continue;
      vertID.irn[k] = pirn[i] + 1;
      vertID.jcn[k] = pjcn[i] + 1;
      if (pirn[i] == pjcn[i]) {
        vertID.a[k]   = pa[i] + 1e-2;
        diagonalVisited[pirn[i]] = true;
      } else
        vertID.a[k]   = pa[i];
      k++;
    }
    for (int i = 0; i < nvar + nconI; i++) {
      if (diagonalVisited[i])
        continue;
      vertID.irn[k] = i + 1;
      vertID.jcn[k] = i + 1;
      vertID.a[k] = 1e-2;
      k++;
      diagonalVisited[i] = true;
    }
    vertID.nz = k;

#ifdef RESTORATIONPRINT
    std::cout << "Augmented Matrix: " << std::endl;
    for (int i = 0; i < vertID.nz; i++) {
      std::cout << '(' << vertID.irn[i] << ','
           << vertID.jcn[i] << ") = "
           << vertID.a[i] << std::endl;
    }
#endif


    Vector dualResidue(*env);
    Vector primalResidue(*env, nvar + nconI);
    pReal rdx = 0;
    pReal rpx = primalResidue.get_doublex();
    Real normrd = 1.0;
    Real normrp = 1.0;

    dualResidue = normGrad;
    rdx = dualResidue.get_doublex();
    //******************************************************
    while ( ( (normrd > 1e-6) || (normrp > 1e-6) || (innerMu > 1e-9) ) && 
            (normc > rho) ) {
      iter++;
      if (iter > maxIterIPRestoration)
        break;
      for (Int i = 0; i < numUpper; i++)
        rdx[upperIndex[i]] += Cx[upperIndex[i]]*multUpper[i];
      for (Int i = 0; i < numLower; i++)
        rdx[lowerIndex[i]] -= Cx[lowerIndex[i]]*multLower[i];
      normrd = dualResidue.norm();

      for (int i = 0; i < nvar + nconI; i++)
        rpx[i] = Cx[i]*dx[i] - slack[i];
      normrp = primalResidue.norm();

      vertID.rhs = new double[nvar + nconI + numUpper + numLower];
      for (int i = 0; i < nvar + nconI; i++)
        vertID.rhs[i] = -rdx[i];
      for (int i = 0; i < numUpper; i++) {
        int j = nvar + nconI + i;
        vertID.rhs[j] = upper[i] - slack[upperIndex[i]] - rpx[upperIndex[i]];
      }
      for (int i = 0; i < numLower; i++) {
        int j = nvar + nconI + numUpper + i;
        vertID.rhs[j] = slack[lowerIndex[i]] - lower[i] + rpx[lowerIndex[i]];
      }

#ifdef RESTORATIONPRINT
      std::cout << "RHS: " << std::endl;
      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        std::cout << vertID.rhs[i] << std::endl;
#endif

      vertID.icntl[0] = -1;
      vertID.icntl[1] = -1;
      vertID.icntl[2] = -1;
      vertID.icntl[3] = 0;
      vertID.job = 6;

      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        rhsClone[i] = vertID.rhs[i];
      //Factorize matrix
      dmumps_c(&vertID);

#ifdef RESTORATIONPRINT
      std::cout << "SOL: " << std::endl;
      for (int i = 0 ; i < nvar + nconI + numUpper + numLower; i++)
        std::cout << vertID.rhs[i] << std::endl;

      Real mumpsResidue = 0.0;
      Real mumpsAxmb[vertID.n];

      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        mumpsAxmb[i] = rhsClone[i];
      for (int k = 0; k < vertID.nz; k++) {
        int i = vertID.irn[k], j = vertID.jcn[k];
        Real Aij = vertID.a[k];
        Real xj = vertID.rhs[j - 1], xi = vertID.rhs[i - 1];
        mumpsAxmb[i - 1] -= Aij*xj;
        if (i != j)
          mumpsAxmb[j - 1] -= Aij*xi;
      }
      for (int k = 0; k < nvar + nconI + numUpper + numLower; k++) {
        mumpsResidue += pow(mumpsAxmb[k], 2);
      }
      std::cout << "||Ax - b||^2 = " << mumpsResidue << std::endl;
#endif


      Real dStep[nvar + nconI], slackStep[nvar + nconI];
      Real upperStep[numUpper], lowerStep[numLower];

      for (int i = 0; i < nvar + nconI; i++) {
        dStep[i] = vertID.rhs[i];
        slackStep[i] = Cx[i]*dStep[i] + rpx[i];
      }
      for (int i = 0; i < numUpper; i++) {
        int j = nvar + nconI + i;
        upperStep[i] = vertID.rhs[j];
      }
      for (int i = 0; i < numLower; i++) {
        int j = nvar + ncon + numUpper + i;
        lowerStep[i] = vertID.rhs[j];
      }

      Real alphaAff = 1.0;
      innerMu = 0.0;
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        innerMu += (upper[i] - slack[j])*multUpper[i];
        if (slackStep[j] > 0)
          alphaAff = Min(alphaAff, (upper[i] - slack[j])/slackStep[j]);
        if (upperStep[i] < 0)
          alphaAff = Min(alphaAff, -multUpper[i]/upperStep[i]);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMu -= (lower[i] - slack[j])*multLower[i];
        if (slackStep[j] < 0)
          alphaAff = Min(alphaAff, (lower[i] - slack[j])/slackStep[j]);
        if (lowerStep[i] < 0)
          alphaAff = Min(alphaAff, -multLower[i]/lowerStep[i]);
      }
      innerMu /= (numUpper + numLower);
      if (alphaAff < 1e-9)
        break;
      
      Real innerMuAff = 0.0;
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        innerMuAff += (upper[i] - slack[j] - alphaAff*slackStep[j])*
                      (multUpper[i] + alphaAff*upperStep[i]);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMuAff -= (lower[i] - slack[j] - alphaAff*slackStep[j])*
                      (multLower[i] + alphaAff*lowerStep[i]);
      }
      innerMuAff /= numUpper + numLower;
      Real sigma = Min(pow(innerMuAff/innerMu, 3), 1e2);

      for (Int i = 0; i < nvar + nconI; i++)
        rdx[i] = normgx[i];
      for (Int i = 0; i < numUpper; i++)
        rdx[upperIndex[i]] += Cx[upperIndex[i]]*multUpper[i];
      for (Int i = 0; i < numLower; i++)
        rdx[lowerIndex[i]] -= Cx[lowerIndex[i]]*multLower[i];
      normrd = dualResidue.norm();

      for (int i = 0; i < nvar + nconI; i++)
        rpx[i] = Cx[i]*dx[i] - slack[i];
      normrp = primalResidue.norm();

      for (int i = 0; i < nvar + nconI; i++)
        vertID.rhs[i] = -rdx[i];
      for (int i = 0; i < numUpper; i++) {
        int j = nvar + nconI + i;
        vertID.rhs[j] = upper[i] - slack[upperIndex[i]] - 
          slackStep[upperIndex[i]]*upperStep[i]/upper[i] -
          sigma*innerMu/multUpper[i] - rpx[upperIndex[i]];
      }
      for (int i = 0; i < numLower; i++) {
        int j = nvar + nconI + numUpper + i;
        vertID.rhs[j] = slack[lowerIndex[i]] - lower[i] +
          slackStep[lowerIndex[i]]*lowerStep[i]/lower[i] -
          sigma*innerMu/multLower[i] + rpx[lowerIndex[i]];
      }

#ifdef RESTORATIONPRINT
      std::cout << "RHS: " << std::endl;
      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        std::cout << vertID.rhs[i] << std::endl;
#endif

      vertID.icntl[0] = -1;
      vertID.icntl[1] = -1;
      vertID.icntl[2] = -1;
      vertID.icntl[3] = 0;
      vertID.job = 6;

      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        rhsClone[i] = vertID.rhs[i];
      dmumps_c(&vertID);

#ifdef RESTORATIONPRINT
      std::cout << "SOL: " << std::endl;
      for (int i = 0 ; i < nvar + nconI + numUpper + numLower; i++)
        std::cout << vertID.rhs[i] << std::endl;

      for (int i = 0; i < nvar + nconI + numUpper + numLower; i++)
        mumpsAxmb[i] = rhsClone[i];
      for (int k = 0; k < vertID.nz; k++) {
        int i = vertID.irn[k], j = vertID.jcn[k];
        Real Aij = vertID.a[k];
        Real xj = vertID.rhs[j - 1], xi = vertID.rhs[i - 1];
        mumpsAxmb[i - 1] -= Aij*xj;
        if (i != j)
          mumpsAxmb[j - 1] -= Aij*xi;
      }
      for (int k = 0; k < nvar + nconI + numUpper + numLower; k++) {
        mumpsResidue += pow(mumpsAxmb[k], 2);
      }
      std::cout << "||Ax - b||^2 = " << mumpsResidue << std::endl;
#endif

      for (int i = 0; i < nvar + nconI; i++) {
        dStep[i] = vertID.rhs[i];
        slackStep[i] = Cx[i]*dStep[i] + rpx[i];
      }
      for (int i = 0; i < numUpper; i++) {
        int j = nvar + nconI + i;
        upperStep[i] = vertID.rhs[j];
      }
      for (int i = 0; i < numLower; i++) {
        int j = nvar + ncon + numUpper + i;
        lowerStep[i] = vertID.rhs[j];
      }

      Real alphaPrimal = 1.0, alphaDual = 1.0;
      Real alphaPrimalMax = 1.0, alphaDualMax = 1.0;
      innerMu = 0.0;
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        innerMu += (upper[i] - slack[upperIndex[i]])*multUpper[i];
        if (slackStep[j] > 0)
          alphaPrimalMax = Min(alphaPrimalMax, (1 - 1e-1)*(upper[i] - slack[j])/slackStep[j]);
        if (upperStep[i] < 0)
          alphaDualMax = Min(alphaDualMax, (1e-1 - 1)*multUpper[i]/upperStep[i]);
        assert(alphaPrimalMax > 0);
        assert(alphaDualMax > 0);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMu -= (lower[i] - slack[lowerIndex[i]])*multLower[i];
        if (slackStep[j] < 0)
          alphaPrimalMax = Min(alphaPrimalMax, (1 - 1e-1)*(lower[i] - slack[j])/slackStep[j]);
        if (lowerStep[i] < 0)
          alphaDualMax = Min(alphaDualMax, (1e-1 - 1)*multLower[i]/lowerStep[i]);
        assert(alphaPrimalMax > 0);
        assert(alphaDualMax > 0);
      }
      innerMu /= (numUpper + numLower);

      /* Armijo test
       * search for t such that
       * f(x + t*d) < f(x) + a*t*dot(g(x),d) 
       * In our situation,
       * f(x) = 0.5*||A*x + h||^2 */
      Real dotgd = 0.0;
      for (int i = 0; i < nvar + nconI; i++)
        dotgd += Cx[i]*dStep[i]*normgx[i];
      Real oldObjFun = 0.5*pow(normc, 2), objFun = oldObjFun;
      alphaPrimal = alphaPrimalMax;
      for (int i = 0; i < nvar; i++)
        xcx[i] = oldxcx[i] + Cx[i]*dx[i] + alphaPrimal*Cx[i]*dStep[i];
      for (int i = 0; i < nconI; i++) {
        int j = nvar + i;
        scx[i] = oldscx[i] + Cx[j]*dx[j] + alphaPrimal*Cx[j]*dStep[j];
      }
      call_ccfsg_xc(dciFalse);
      normc = c->norm();
      objFun = 0.5*pow(normc, 2);

      Vector oldDualResidue(dualResidue), oldPrimalResidue(primalResidue);
      pReal oldrdx = oldDualResidue.get_doublex();
      pReal oldrpx = oldPrimalResidue.get_doublex();

      for (int i = 0; i < nvar + nconI; i++) {
        rpx[i] = oldrpx[i] + alphaPrimal*(Cx[i]*dStep[i] - slackStep[i]);
      }

      Vector dStepVector(*env, nvar + nconI, dStep);
      dStepVector.sdmult(*J, 0, one, zero, dStepVector);
      dStepVector.sdmult(*J, 1, one, zero, dStepVector);
      pReal dsvx = dStepVector.get_doublex();

      for (int i = 0; i < nvar + nconI; i++) {
        rdx[i] = oldrdx[i] + alphaPrimal*dsvx[i];
      }
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        rdx[j] += Min(alphaPrimal, alphaDualMax)*Cx[j]*upperStep[i];
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        rdx[j] -= Min(alphaPrimal, alphaDualMax)*Cx[j]*lowerStep[i];
      }

      //IntervalMinimization
      Real alphaPrimalLeft = 0, alphaPrimalRight = alphaPrimalMax;
      Real leftObjFun, rightObjFun;
      leftObjFun = InteriorPointObjFun (alphaPrimalLeft, 
          Min(alphaDualMax, alphaPrimalLeft), numUpper, numLower, oldrdx, oldrpx,
          dsvx, Cx, dStep, slackStep, upperStep, lowerStep, upperIndex,
          lowerIndex);
      rightObjFun = InteriorPointObjFun (alphaPrimalRight, 
          Min(alphaDualMax,alphaPrimalRight), numUpper, numLower, oldrdx, oldrpx,
          dsvx, Cx, dStep, slackStep, upperStep, lowerStep, upperIndex,
          lowerIndex);
      while (rightObjFun - leftObjFun > alphaPrimalMax/2e8) {
        Real x = (alphaPrimalLeft + alphaPrimalRight)/2;
        Real fx = InteriorPointObjFun (x, 
          Min(alphaDualMax, x), numUpper, numLower, oldrdx, oldrpx,
          dsvx, Cx, dStep, slackStep, upperStep, lowerStep, upperIndex,
          lowerIndex);

        if (leftObjFun < rightObjFun) {
          alphaPrimalRight = x;
          rightObjFun = fx;
        } else {
          alphaPrimalLeft = x;
          leftObjFun = fx;
        }
      }
      alphaPrimal = (alphaPrimalLeft + alphaPrimalRight)/2;
      for (int i = 0; i < nvar; i++)
        xcx[i] = oldxcx[i] + Cx[i]*dx[i] + alphaPrimal*Cx[i]*dStep[i];
      for (int i = 0; i < nconI; i++) {
        int j = nvar + i;
        scx[i] = oldscx[i] + Cx[j]*dx[j] + alphaPrimal*Cx[j]*dStep[j];
      }
      //IntervalMinimization End
/*       while (objFun > oldObjFun) {
 *         alphaPrimal /= 2;
 *         if (alphaPrimal < alphaPrimalMax*2e-8)
 *           break;
 *         for (int i = 0; i < nvar; i++)
 *           xcx[i] = oldxcx[i] + Cx[i]*dx[i] + alphaPrimal*Cx[i]*dStep[i];
 *         for (int i = 0; i < nconI; i++) {
 *           int j = nvar + i;
 *           scx[i] = oldscx[i] + Cx[j]*dx[j] + alphaPrimal*Cx[j]*dStep[j];
 *         }
 *         call_ccfsg_xc(dciFalse);
 *         normc = c->norm();
 *         objFun = 0.5*pow(normc, 2);
 *       }
 */
      for (int i = 0; i < nvar; i++) {
        dx[i] += alphaPrimal*dStep[i];
        slack[i] += alphaPrimal*slackStep[i];
      }
      for (int i = 0; i < nconI; i++) {
        int j = nvar + i;
        dx[j] += alphaPrimal*dStep[j];
        slack[j] += alphaPrimal*slackStep[j];
      }
      alphaDual = Min(alphaDualMax, alphaPrimal);

      //This really belongs here?
/*       numLower = 0;
 *       numUpper = 0;
 *       for (Int i = 0; i < nvar; i++) {
 *         Real zi = xcx[i], li = blx[i], ui = bux[i];
 *         if (li > -dciInf) {
 *           lower[numLower] = (li - zi) * (1 - 1e-1);
 *           lowerIndex[numLower++] = i;
 *         }
 *         if (ui < dciInf) {
 *           upper[numUpper] = (ui - zi) * (1 - 1e-1);
 *           upperIndex[numUpper++] = i;
 *         }
 *       }
 *       for (Int i = 0; i < nconI; i++) {
 *         Int j = nvar + i;
 *         Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
 *         if (li > -dciInf) {
 *           lower[numLower] = (li - zi) * (1 - 1e-1);
 *           lowerIndex[numLower++] = j;
 *         }
 *         if (ui < dciInf) {
 *           upper[numUpper] = (ui - zi) * (1 - 1e-1);
 *           upperIndex[numUpper++] = j;
 *         }
 *       }
 */
      //^^
#ifdef RESTORATIONPRINT
      std::cout << "xc = " << std::endl;
      xc->print_more();
      std::cout << "sc = " << std::endl;
      sc->print_more();
      std::cout << "|c| = " << c->norm() << std::endl;
#endif

      checkInfactibility ();
      for (int i = 0; i < numUpper; i++) {
        multUpper[i] += alphaDual*upperStep[i];
        assert(multUpper[i] > 0);
        assert(slack[upperIndex[i]] < upper[i]);
      }
      for (int i = 0; i < numLower; i++) {
        multLower[i] += alphaDual*lowerStep[i];
        assert(multLower[i] > 0);
        assert(slack[lowerIndex[i]] > lower[i]);
      }

      call_ccfsg_xc(dciFalse);
      normc = c->norm();

      if ( ( (normrp < 1e-6) && (normrd < 1e-6) && (innerMu < 1e-6) ) ||
           (normc < rho) )
        break;
      
      for (int i = 0; i < numUpper; i++) {
        vertID.irn[2*i] = nvar + nconI + i + 1;
        vertID.jcn[2*i] = i + 1;
        vertID.a[2*i]   = Cx[upperIndex[i]];

        vertID.irn[2*i + 1] = nvar + ncon + i + 1;
        vertID.jcn[2*i + 1] = nvar + ncon + i + 1;
        vertID.a[2*i + 1]   = -(upper[i] - slack[upperIndex[i]])/multUpper[i];
      }
      for (int i = 0; i < numLower; i++) {
        vertID.irn[2*i + 2*numUpper] = nvar + nconI + numUpper + i + 1;
        vertID.jcn[2*i + 2*numUpper] = i + 1;
        vertID.a[2*i + 2*numUpper]   = -Cx[lowerIndex[i]];

        vertID.irn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
        vertID.jcn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
        vertID.a[2*i + 1 + 2*numUpper]   = (lower[i] - slack[lowerIndex[i]])/multLower[i];
      }

#ifdef RESTORATIONPRINT
      std::cout << "Augmented Matrix: " << std::endl;
      for (int i = 0; i < vertID.nz; i++) {
        std::cout << '(' << vertID.irn[i] << ','
             << vertID.jcn[i] << ") = "
             << vertID.a[i] << std::endl;
      }
#endif

      normGrad = *c;
      normGrad.sdmult(*J, 0, one, one, d);
      normGrad.sdmult(*J, 1, one, zero, normGrad);
      normgx = normGrad.get_doublex();

      dualResidue = normGrad;
      rdx = dualResidue.get_doublex();
    }

    if (normc > oldnormc) {
      *xc = oldxc;
      *sc = oldsc;
      normc = oldnormc;
      DeltaV /= 4;
      call_ccfsg(dciFalse);
    }

    free(vertID.rhs);
    free(vertID.irn);
    free(vertID.jcn);
    free(vertID.a);

    free(slack);
    free(multUpper);
    free(multLower);

    return 0;
  }
}
