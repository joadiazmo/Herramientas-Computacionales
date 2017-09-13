#include<iostream>
#include<cmath>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>

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
  evals = eigensolver.eigenvalues().real();
  eivec = eigensolver.eigenvectors();
  
  std::sort(evals.data(), evals.data() + evals.size());
  
  return evals(index);
}





#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>

void set_as_hilbert(Eigen::MatrixXd & M);
double condition_number(Eigen::MatrixXd & M);

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  
  for (int N = 1; N <= 4096; N *= 2)
    {
      Eigen::MatrixXd M(N, N);
      set_as_hilbert(M);
      std::cout << N << " " << condition_number(M) << std::endl;
    }
  return 0;
}

void set_as_hilbert(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n)
    {
      for (int m = 0; m < M.cols(); ++m)
	{
	  M(n, m) = 1./(n+m+1.);
	}
    }
}

double condition_number(Eigen::MatrixXd & M)
{
  // Escriba acá el código para calcular el condition number. Puede
  // ser una sola linea! Se le recomienda buscar como calcular la
  // norma de una matriz de eigen (funcion .norm() ?)
}
