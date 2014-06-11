#include <iostream>
#include "dci.h"
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace DCI;

Real rand_between (Real a, Real b) {
  Real x = (rand()%1001)/1000.0;
  return (b - a)*x + a;
}

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Int npt = (*n)/2;
  Real sum;

  *f = 0.0;

  if (*grad == dciTrue) {
    for (Int i = 0; i < *n; i++)
      g[i] = 0.0;
  }

  for (Int i = 0; i < npt - 1; i++) {
    for (Int j = i+1; j < npt; j++) {
      sum = pow(x[2*i] - x[2*j], 2) +
            pow(x[2*i+1] - x[2*j+1], 2);
      *f += 1/sum;
      sum *= sum;

      if (*grad == dciTrue) {
        g[2*i] -= 2 * ( x[2*i] - x[2*j] ) / sum;
        g[2*j] += 2 * ( x[2*i] - x[2*j] ) / sum;
        g[2*i + 1] -= 2 * ( x[2*i + 1] - x[2*j + 1] ) / sum;
        g[2*j + 1] += 2 * ( x[2*i + 1] - x[2*j + 1] ) / sum;
      }
    }
  }
}

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Int npt = (*n)/2;
  Real sum, sum2, sum4;
  Int nnz = (*n)*(*n + 1)/2;

  Int I[nnz], J[nnz];
  Real H[nnz];

  Int k = 0;
  for (Int i = 0; i < *n; i++) {
    for (Int j = 0; j <= i; j++) {
      I[k] = i;
      J[k] = j;
      k++;
    }
  }

  for (Int i = 0; i < nnz; i++)
    H[i] = 0.0;

  for (Int i = 0; i < *n; i++)
    q[i] = 0.0;

  for (Int i = 0; i < npt - 1; i++) {
    for (Int j = i+1; j < npt; j++) {
      sum = pow(x[2*i] - x[2*j], 2) + pow(x[2*i+1] - x[2*j+1], 2);
      sum2 = sum*sum;
      sum4 = sum2*sum2;

      Int a = 2*i, b = 2*i+1, c = 2*j, d = 2*j+1;

      H[d*(d-1)/2 + d] += (-2 * sum2 + 8 * (x[b] - x[d]) * sum * (x[b] - x[d]) )/sum4;
      H[d*(d-1)/2 + c] +=            ( 8 * (x[b] - x[d]) * sum * (x[a] - x[c]) )/sum4;
      H[d*(d-1)/2 + b] += ( 2 * sum2 - 8 * (x[b] - x[d]) * sum * (x[b] - x[d]) )/sum4;
      H[d*(d-1)/2 + a] +=            (-8 * (x[b] - x[d]) * sum * (x[a] - x[c]) )/sum4;
      H[c*(c-1)/2 + c] += (-2 * sum2 + 8 * (x[a] - x[c]) * sum * (x[a] - x[c]) )/sum4;
      H[c*(c-1)/2 + b] +=            (-8 * (x[a] - x[c]) * sum * (x[b] - x[d]) )/sum4;
      H[c*(c-1)/2 + a] += ( 2 * sum2 - 8 * (x[a] - x[c]) * sum * (x[a] - x[c]) )/sum4;
      H[b*(b-1)/2 + b] += (-2 * sum2 + 8 * (x[b] - x[d]) * sum * (x[b] - x[d]) )/sum4;
      H[b*(b-1)/2 + a] +=            ( 8 * (x[b] - x[d]) * sum * (x[a] - x[c]) )/sum4;
      H[a*(a-1)/2 + a] += (-2 * sum2 + 8 * (x[a] - x[c]) * sum * (x[a] - x[c]) )/sum4;
    }
  }

  for (k = 0; k < nnz; k++) {
    Int i = I[k], j = J[k], hij = H[k];
    q[i] += hij * p[j];
    q[j] += hij * p[i];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Int npt = (*n)/2;
  Real sum;

  *f = 0.0;

  for (Int i = 0; i < npt; i++) {
    for (Int j = i+1; j < npt; j++) {
      sum = 1/( pow(x[2*i] - x[2*j], 2) +
                pow(x[2*i+1] - x[2*j+1], 2) );
      *f += sum;
    }
    c[2*i] = x[2*i+1] - x[2*i] - 0.5;
    c[2*i+1] = x[2*i+1] - x[2*i] + 0.5;
  }
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Int npt = (*n)/2;
  Real sum;

  for (Int i = 0; i < npt; i++) {
    c[2*i] = x[2*i+1] - x[2*i] - 0.5;
    c[2*i+1] = x[2*i+1] - x[2*i] + 0.5;
  }
  if (*Grad == dciFalse)
    return;

  for (Int i = 0; i < npt; i++) {
    indvar[4*i] = 2*i;
    indfun[4*i] = 2*i;
    indvar[4*i + 1] = 2*i;
    indfun[4*i + 1] = 2*i + 1;
    indvar[4*i + 2] = 2*i + 1;
    indfun[4*i + 2] = 2*i;
    indvar[4*i + 3] = 2*i + 1;
    indfun[4*i + 3] = 2*i + 1;
    J[4*i] = -1;
    J[4*i + 1] = -1;
    J[4*i + 2] =  1;
    J[4*i + 3] =  1;
  }
  *nnzJ = 4*npt;
}

int main () {
  Int npt = 100;
  Int n = 2*npt, m = 2*npt;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  srand(time(0));

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);
//  dci.set_ccifg (CCIFG);

  for (Int i = 0; i < n; i++) {
    x[i] = rand_between(-1,1);
    bl[i] = -1;
    bu[i] = 1;
  }

  for (Int i = 0; i < npt; i++) {
    equatn[2*i] = dciFalse;
    equatn[2*i+1] = dciFalse;
    y[2*i] = 0;
    y[2*i+1] = 0;
    // For infeasible points, the 4 next lines must be uncommented
    cl[2*i] = 0;
    cu[2*i] = dciInf;
    cl[2*i+1] = -dciInf;
    cu[2*i+1] = 0;
    // For feasible points, the 4 next lines must be uncommented
//    cl[2*i+1] = 0;
//    cu[2*i+1] = dciInf;
//    cl[2*i] = -dciInf;
//    cu[2*i] = 0;
  }

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_equatn (m, equatn);

  dci.start ();
  dci.solve ();
//  dci.show();

  pReal sol = dci.get_x();
  std::cout << "x = [";
  for (Int i = 0; i < npt; i++) {
    std::cout << sol[2*i] << ' ';
  }
  std::cout << "];\n";
  std::cout << "y = [";
  for (Int i = 0; i < npt; i++) {
    std::cout << sol[2*i+1] << ' ';
  }
  std::cout << "];\n";
  std::cout << "figure" << std::endl;
  std::cout << "plot(x,y,'ro')" << std::endl;
  std::cout << "hold on; fplot(@(t) t + 0.5,[-1 1]); fplot(@(t) t - 0.5,[-1,1]);" << std::endl;
  std::cout << "axis([-1 1 -1 1])" << std::endl;

}
