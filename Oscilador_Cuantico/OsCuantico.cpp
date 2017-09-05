#include<iostream>
#include<cmath>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>
//#include<eigen3/Eigen/Dense>
//#include<eigen3/Eigen/EigenValues>

using namespace Eigen;

void set_H0(Eigen::MatrixXd & M);
void set_X(Eigen::MatrixXd & M);
double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, const double lambda, const int index);

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  
  double lambda = 0.2;
  
  for (int N = 2; N <= 1024; N *= 2)
    {
      Eigen::MatrixXd X(N, N), H0(N, N), H(N, N);
      set_H0(H0);
      set_X(X);
      std::cout << 1.0/N << " " << eigen_energy(H0, X, lambda, 0) << std::endl;
    }
  
  return 0;
}

void set_H0(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n)
    {
      M(n, n) += (n+1/2.0);
    }
}

void set_X(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n)
    {
    for (int m = 0; m < M.cols(); ++m)
      {
	if (n == m+1) M(n, m) += std::sqrt((m+1.0)/2.0);
	if (n == m-1) M(n, m) += std::sqrt((m)/2.0);
      }
    }
}

double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, const double lambda, const int index)
{
  Eigen::MatrixXd A, eivec;
  Eigen::VectorXd evals;
  
  A = H + lambda*X*X*X*X;
  
  SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
  
  if (eigensolver.info() != Success) abort();
  evals = eigensolver.eigenvalues();
  eivec = eigensolver.eigenvectors();
  
  std::sort(evals.data(), evals.data() + evals.size());
  
  return evals(index);
}
