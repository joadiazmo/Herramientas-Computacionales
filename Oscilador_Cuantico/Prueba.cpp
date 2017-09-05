#include<iostream>
#include<cmath>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>
//#include<eigen3/Eigen/EigenValues>

using Eigen::MatrixXd;

void set_H0(Eigen::MatrixXd & M);
int N = 3;

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  
  double lambda = 0.2;
  
  Eigen::MatrixXd H0(N, N);
  set_H0(H0);
  
  std::cout << H0  << std::endl;
  
  return 0;
}

void set_H0(Eigen::MatrixXd & A)
{
  A.setZero();
  for (int n = 0; n < A.cols(); ++n)
    {
      A(n, n) = (n+1/2.0);
    }
}
