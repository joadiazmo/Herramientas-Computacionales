#include<iostream>
#include<cmath>
#include<algorithm>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Core>

using namespace Eigen;

void set_as_hilbert(MatrixXd & M);
double condition_number(Eigen::MatrixXd & M);

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  for (int N = 1; N <= 4096; N *= 2)
    {
      MatrixXd M(N, N);
      set_as_hilbert(M);
      std::cout << N << " " << condition_number(M) << std::endl;
    }
  return 0;
}

void set_as_hilbert(MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n)
    {
      for (int m = 0; m < M.cols(); ++m)
	{
	  M(n, m) += 1./(m+n+1.);
	}
    }
}


double condition_number(MatrixXd & M)
{
  double k;
  MatrixXd I;
  
  I = M.inverse();
  k = M.norm() * I.norm();
  return k;
}
