#include<iostream>
#include<cmath>
#include<algorithm>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Core>

using namespace Eigen;

void set_H(MatrixXd & M);
double norma;
double determinante;

int main(int argc, char **argv)
{
  Matrix2d H, I;

  H << 1, 2,
    3, 4;
  
  I = H.inverse();
  
  norma = H.norm();
  determinante = H.determinant();
  
  std::cout << H << "\n\n"  << I  << "\n\n" << norma << "\t"  << determinante  << std::endl;
  
  return 0;
}
