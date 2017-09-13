#include<iostream>
#include<cmath>
#include<armadillo>
//#include<algorithm>

using namespace std;
using namespace arma;

void set_H0(mat & M);
void set_X(mat & M);
double eigen_energy(mat & H, mat & X, const double lambda, const int index);


int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  
  double lambda = 0.2;
  
  for (int N = 2; N <= 1024; N *= 2)
    {
      mat X(N, N), H0(N, N), H(N, N);
      set_H0(H0);
      set_X(X);
      cout << 1.0/N << " " << eigen_energy(H0, X, lambda, 0) <<endl;
    }
  
  return 0;
}

void set_H0(mat & M)
{
  M.fill(0.0);
  for (int n = 0; n < M.n_cols; ++n)
    {
      M(n, n) += (n+1/2.0);
    }
}

void set_X(mat & M)
{
  M.fill(0.0);
  for (int n = 0; n < M.n_cols; ++n)
    {
    for (int m = 0; m < M.n_cols; ++m)
      {
	if (n == m+1) M(n, m) += std::sqrt((m+1.0)/2.0);
	if (n == m-1) M(n, m) += std::sqrt((m)/2.0);
      }
    }
}

double eigen_energy(mat & H, mat & X, const double lambda, const int index)
{
  mat A, eivec;
  vec evals;
  
  A = H + lambda*X*X*X*X;
    
  eigs_sym(evals, eivec, A, 2);
  
  return evals(index);
}
