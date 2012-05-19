#include "interface.h"
#include <cmath>

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
  Int Interface::InteriorPointRestoration () {
    Vector d(*env);
    d.reset(nvar + nconI, 0.0);
    Real *slack, *multUpper, *multLower;
    Real upper[nvar + nconI], lower[nvar + nconI]; //Bounds for C*d.
    Vector matrixC(*env);
    matrixC.reset(nvar + nconI, 1.0);
    Int iter = 0, maxIterIPRestoration = 100;

    Int numUpper = 0, numLower = 0;
    Int upperIndex[nvar + nconI], lowerIndex[nvar + nconI];

    if (ScaleVertical)
      scale_xc(matrixC);

    pReal Cx = matrixC.get_doublex();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], li = blx[i], ui = bux[i];
      if (li > -dciInf) {
        lower[numLower] = (li - zi) * (1 - epsmu);
        lowerIndex[numLower++] = i;
      }
      if (ui < dciInf) {
        upper[numUpper] = (ui - zi) * (1 - epsmu);
        upperIndex[numUpper++] = i;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      if (li > -dciInf) {
        lower[numLower] = (li - zi) * (1 - epsmu);
        lowerIndex[numLower++] = j;
      }
      if (ui < dciInf) {
        upper[numUpper] = (ui - zi) * (1 - epsmu);
        upperIndex[numUpper++] = j;
      }
    }

    Real innerMu = 1.0;
    multUpper = new Real[numUpper];
    multLower = new Real[numLower];
    slack     = new Real[nvar + nconI];

    for (int i = 0; i < numUpper; i++)
      multUpper[i] = 1.0;
    for (int i = 0; i < numLower; i++)
      multLower[i] = 1.0;
    for (int i = 0; i < nvar + nconI; i++)
      slack[i] = 0.0;

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
    for (int i = 0; i < gtnnz; i++) {
      if (pirn[i] < pjcn[i])
        continue;
      vertID.irn[k] = pirn[i] + 1;
      vertID.jcn[k] = pjcn[i] + 1;
      vertID.a[k]   = pa[i];
      k++;
    }
    vertID.nz = k;

    std::cout << "Augmented Matrix: " << std::endl;
    for (int i = 0; i < vertID.nz; i++) {
      std::cout << '(' << vertID.irn[i] << ','
           << vertID.jcn[i] << ") = "
           << vertID.a[i] << std::endl;
    }
    Vector dualResidue(*env);
    Vector normGrad(*env);
    Vector primalResidue(*env, nvar + nconI);
    pReal rdx = 0;
    pReal rpx = primalResidue.get_doublex(), dx = d.get_doublex();
    pReal normgx = 0;
    Real one[2] = {1,0}, zero[2] = {0,0};
    Real normrd = 1.0;
    Real normrp = 1.0;

    normGrad.sdmult(*J, 1, one, zero, *c);
    normgx = normGrad.get_doublex();
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

      vertID.icntl[0] = -1;
      vertID.icntl[1] = -1;
      vertID.icntl[2] = -1;
      vertID.icntl[3] = 0;
      vertID.job = 6;

      //Factorize matrix
      dmumps_c(&vertID);

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
      Real sigma = pow(innerMuAff/innerMu, 3);

      for (Int i = 0; i < nvar + nconI; i++)
        rdx[i] = normgx[i];
      for (Int i = 0; i < numUpper; i++)
        rdx[upperIndex[i]] += Cx[upperIndex[i]]*multUpper[i];
      for (Int i = 0; i < numLower; i++)
        rdx[lowerIndex[i]] -= Cx[lowerIndex[i]]*multLower[i];

      for (int i = 0; i < nvar + nconI; i++)
        rpx[i] = Cx[i]*dx[i] - slack[i];

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

      vertID.icntl[0] = -1;
      vertID.icntl[1] = -1;
      vertID.icntl[2] = -1;
      vertID.icntl[3] = 0;
      vertID.job = 6;

      dmumps_c(&vertID);

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
      innerMu = 0.0;
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        innerMu += (upper[i] - slack[upperIndex[i]])*multUpper[i];
        if (slackStep[j] > 0)
          alphaPrimal = Min(alphaPrimal, (1 - epsmu)*(upper[i] - slack[j])/slackStep[j]);
        if (upperStep[i] < 0)
          alphaDual = Min(alphaDual, (epsmu - 1)*multUpper[i]/upperStep[i]);
        assert(alphaPrimal > 0);
        assert(alphaDual > 0);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMu -= (lower[i] - slack[lowerIndex[i]])*multLower[i];
        if (slackStep[j] < 0)
          alphaPrimal = Min(alphaPrimal, (1 - epsmu)*(lower[i] - slack[j])/slackStep[j]);
        if (lowerStep[i] < 0)
          alphaDual = Min(alphaDual, (epsmu - 1)*multLower[i]/lowerStep[i]);
        assert(alphaPrimal > 0);
        assert(alphaDual > 0);
      }
      innerMu /= (numUpper + numLower);

      for (int i = 0; i < nvar; i++) {
        dx[i] += alphaPrimal*dStep[i];
        slack[i] += alphaPrimal*slackStep[i];
        xcx[i] += alphaPrimal*Cx[i]*dStep[i];
      }
      for (int i = 0; i < nconI; i++) {
        int j = nvar + i;
        dx[j] += alphaPrimal*dStep[j];
        slack[j] += alphaPrimal*slackStep[j];
        scx[i] += alphaPrimal*Cx[j]*dStep[j];
      }
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

      call_ccfsg(dciFalse);
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

      normGrad = *c;
      normGrad.sdmult(*J, 0, one, one, d);
      normGrad.sdmult(*J, 1, one, zero, normGrad);
      normgx = normGrad.get_doublex();

      dualResidue = normGrad;
      rdx = dualResidue.get_doublex();
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
