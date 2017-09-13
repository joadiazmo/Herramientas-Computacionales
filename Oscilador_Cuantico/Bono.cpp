#include <iostream>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

using namespace Eigen;

const int N = 1024;
double l = 0.0;

void set_H0(Eigen::MatrixXd & M);
void set_X(Eigen::MatrixXd & M);
double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, double lambda, const int index);

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);

  Eigen::MatrixXd X(N, N), H0(N, N);
  set_H0(H0);
  set_X(X);
  
  for (int i = 0; i <= 10; i++)
    {
      l = 0.1*i;
      
      std::cout << l  << " " << eigen_energy(H0, X, l, 0) << std::endl;
    }
  return 0;
}

void set_H0(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n){
    for (int m = 0; m < M.cols(); ++m){
      if (n == m) M(n, m) += (n+(1/2.));
      
    }
  }
  // Escriba aca el codigo que crea la matriz H0
  // Puede guiarse por la funcion set_X
}

void set_X(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n){
    for (int m = 0; m < M.cols(); ++m){
      if (n == m+1) M(n, m) += std::sqrt((m+1.0)/2.0);
      if (n == m-1) M(n, m) += std::sqrt((m)/2.0);
    }
  }
}

double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, double lambda, const int index)
{
  
  Eigen::MatrixXd Hlambda;
  Hlambda = H + lambda*X*X*X*X;
  
  Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(Hlambda);
  if (eigensolver.info() != Eigen::Success) abort();
  // Implemente aca el calculo de los valores propios, usando la libreria eigen
  // - Calculo de Hlambda :
  // - Calculo de los valores propios (y vectores propios) :
  
  Eigen::VectorXd evals;
  evals = eigensolver.eigenvalues().real();
  // - Extraer los valores propios al vector Eigen::VectorXd evals, solamente
  // la parte real (cuando pida los valores propios escriba .eigenvalues.real()) :
  // ordenar los valores propios :
  std::sort(evals.data(), evals.data() + evals.size());
  // retornar el valor propio (para el estado base index == 0) :
  return evals(index);
}
