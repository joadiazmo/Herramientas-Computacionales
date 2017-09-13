#include<iostream>
//#include<cmath>
//#include<algorithm>
#include<armadillo>

using namespace std;
using namespace arma;

double norma;
double determinante;

int main()
{
  mat H;
  mat I;

  H << 1 << 2 << endr
    << 3 << 4 << endr;
  
  I = inv(H);

  cx_vec eigval = eig_gen(H);
  
  norma = norm(H,2);
  determinante = det(H);
  
  cout << H << "\n\n"  << I  << "\n\n" << norma << "\t"  << determinante  <<endl;
  
  return 0;
}

