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
    Vector d(*env, nvar + nconI, 0.0);
    Real *slack, *multUpper, *multLower;
    Real upper[nvar + nconI], lower[nvar + nconI]; //Bounds for C*d.
    Vector matrixC(*env);
    matrixC.reset(nvar + nconI, 1.0);

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
    for (Int i = 0; i < nvar; i++) {
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
    multUpper = new Real[numUpper](1.0);
    multLower = new Real[numUpper](1.0);
    slack     = new Real[nvar + nconI](0.0);

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
    G.aat (G); // G = G*G' = A'*A;
    Triplet GT(*env);
    T.sparse_to_triplet(G);
    //Create initial matrix
    vertID.nz  = (int)(G.nnz()) + 2*(numUpper + numLower);

    vertID.irn = new int[vertID.nz];
    vertID.jcn = new int[vertID.nz];
    vertID.a   = new double[vertID.nz];

    for (int i = 0; i < numUpper; i++) {
      vertID.irn[2*i] = nvar + nconI + i + 1;
      vertID.jcn[2*i] = i + 1;
      vertiD.a[2*i]   = matrixC[upperIndex[i]];

      vertID.irn[2*i + 1] = nvar + ncon + i + 1;
      vertID.jcn[2*i + 1] = nvar + ncon + i + 1;
      vertID.a[2*i + 1]   = -(upper[i] - slack[upperIndex[i]])/multUpper[i];
    }
    for (int i = 0; i < numLower; i++) {
      vertID.irn[2*i + 2*numUpper] = nvar + nconI + numUpper + i + 1;
      vertID.jcn[2*i + 2*numUpper] = i + 1;
      vertiD.a[2*i + 2*numUpper]   = -matrixC[lowerIndex[i]];

      vertID.irn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
      vertID.jcn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
      vertID.a[2*i + 1 + 2*numUpper]   = (lower[i] - slack[lower[i]])/multLower[i];
    }
    int k = 2*numUpper + 2*numLower;
    int gtnnz = (*GT.get_pnnz());
    int *pirn = GT.get_linti(), *pjcn = GT.get_lintj();
    double *pa = GT.get_doublex();
    for (int i = 0; i < gtnnz; i++) {
      vertID.irn[k + i] = pirn[i] + 1;
      vertID.jcn[k + i] = pjcn[i] + 1;
      vertID.a[k + i]   = pa[i];
    }
    vertID.nz = k + gtnnz;
    Vector dualResidue(*env);
    pReal rdx = dualResidue.get_doublex();
    Vector primalResidue(*env, nvar + nconI);
    pReal rpx = primalResidue.get_doublex(), pReal dx = d.get_doublex();
    Real one[2] = {1,0}, zero[2] = {0,0};

    dualResidue.sdmult(*J, 1, one, zero, *c);
    //******************************************************
    while ( (innerMu > 1e-9) && (normc > rho) ) {
      for (Int i = 0; i < numUpper; i++)
        rdx[i] += Cx[upperIndex[i]]*multUpper[i];
      for (Int i = 0; i < numLower; i++)
        rdx[i] -= Cx[lowerIndex[i]]*multLower[i];

      for (int i = 0; i < nvar + nconI; i++)
        rpx[i] = Cx[i]*dx[i] - slack[i];

      vertID.rhs = new int[nvar + nconI + numUpper + numLower];
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
        innerMu += (upper[i] - slack[upperIndex[i]])*multUpper[i];
        if (slackStep[j] > 0)
          alphaAff = Min(alphaAff, (upper[i] - slack[j])/slackStep[j]);
        if (upperStep[i] > 0)
          alphaAff = Min(alphaAff, -multUpper[i]/upperStep[i]);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMu += (lower[i] - slack[lower[i]])*multLower[i];
        if (slackStep[j] < 0)
          alphaAff = Min(alphaAff, (lower[i] - slack[j])/slackStep[j]);
        if (lowerStep[i] > 0)
          alphaAff = Min(alphaAff, -multLower[i]/lowerStep[i]);
      }
      innerMu /= (numUpper + numLower);
      
      Real innerMuAff = 0.0;
      for (int i = 0; i < numUpper; i++) {
        int j = upperIndex[i];
        inneqMuAff += (upper[i] - slack[j] - alphaAff*slackStep[j])*
                      (multUpper[i] + alphaAff*upperStep[i]);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        inneqMuAff += (lower[i] - slack[j] - alphaAff*slackStep[j])*
                      (multLower[i] + alphaAff*lowerStep[i]);
      }
      Real sigma = pow(innerMuAff/innerMu, 3);

      for (Int i = 0; i < numUpper; i++)
        rdx[i] += Cx[upperIndex[i]]*multUpper[i];
      for (Int i = 0; i < numLower; i++)
        rdx[i] -= Cx[lowerIndex[i]]*multLower[i];

      for (int i = 0; i < nvar + nconI; i++)
        rpx[i] = Cx[i]*dx[i] - slack[i];

      vertID.rhs = new int[nvar + nconI + numUpper + numLower];
      for (int i = 0; i < nvar + nconI; i++)
        vertID.rhs[i] = -rdx[i];
      for (int i = 0; i < numUpper; i++) {
        int j = nvar + nconI + i;
        vertID.rhs[j] = upper[i] - slack[upperIndex[i]] - 
          slackStep[upperIndex[i]]*upper[i] -
          sigma*vertMu/multUpper[i] - rpx[upperIndex[i]];
      }
      for (int i = 0; i < numLower; i++) {
        int j = nvar + nconI + numUpper + i;
        vertID.rhs[j] = slack[lowerIndex[i]] - lower[i] +
          slackStep[lowerIndex[i]]*lower[i] -
          sigma*vertMu/multLower[i] + rpx[lowerIndex[i]];
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
        if (upperStep[i] > 0)
          alphaDual = Min(alphaDual, (epsmu - 1)*multUpper[i]/upperStep[i]);
      }
      for (int i = 0; i < numLower; i++) {
        int j = lowerIndex[i];
        innerMu += (lower[i] - slack[lower[i]])*multLower[i];
        if (slackStep[j] < 0)
          alphaPrimal = Min(alphaPrimal, (1 - epsmu)*(lower[i] - slack[j])/slackStep[j]);
        if (lowerStep[i] > 0)
          alphaDual = Min(alphaDual, (epsmu - 1)*multLower[i]/lowerStep[i]);
      }
      innerMu /= (numUpper + numLower);

      for (int i = 0; i < nvar; i++) {
        dx[i] += alphaPrimal*dStep[i];
        slack[i] += alphaPrimal*slackStep[i];
        xcx[i] += alphaPrimal*dStep[i];
      }
      for (int i = 0; i < nconI; i++) {
        int j = nvar + i;
        dx[j] += alphaPrimal*dStep[j];
        slack[j] += alphaPrimal*slackStep[j];
        scx[i] += alphaPrimal*dStep[j];
      }
      for (int i = 0; i < numUpper; i++) {
        multUpper[i] += alphaDual*upperStep[i];
      }
      for (int i = 0; i < numLower; i++) {
        multLower[i] += alphaDual*lowerStep[i];
      }

      call_ccfsg(dciFalse);
      normc = c->norm();

      if ( (normc < rho) || (innerMu < 1e-6) )
        break;
      
      for (int i = 0; i < numUpper; i++) {
        vertID.irn[2*i] = nvar + nconI + i + 1;
        vertID.jcn[2*i] = i + 1;
        vertiD.a[2*i]   = matrixC[upperIndex[i]];

        vertID.irn[2*i + 1] = nvar + ncon + i + 1;
        vertID.jcn[2*i + 1] = nvar + ncon + i + 1;
        vertID.a[2*i + 1]   = -(upper[i] - slack[upperIndex[i]])/multUpper[i];
      }
      for (int i = 0; i < numLower; i++) {
        vertID.irn[2*i + 2*numUpper] = nvar + nconI + numUpper + i + 1;
        vertID.jcn[2*i + 2*numUpper] = i + 1;
        vertiD.a[2*i + 2*numUpper]   = -matrixC[lowerIndex[i]];

        vertID.irn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
        vertID.jcn[2*i + 1 + 2*numUpper] = nvar + ncon + numUpper + i + 1;
        vertID.a[2*i + 1 + 2*numUpper]   = (lower[i] - slack[lower[i]])/multLower[i];
      }

      dualResidue = *c;
      dualResidue.sdmult(*J, 0, one, one, d);
      dualResidue.sdmult(*J, 1, one, zero, dualResidue);

    }

    free(vertID.rhs);
    free(vertID.irn);
    free(vertID.jcn);
    free(vertID.a);

    free(slack);
    free(multUpper);
    free(multLower);
  }
}
