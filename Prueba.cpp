#include <iostream>
#include <cmath>

using namespace std;

//t es el numero de terminos

int t = 10;
double suma = 0;
double termino = 0;

double f(int n);

int main()
{
  std::cout.precision(16);  std::cout.setf(std::ios::scientific);
  
  for(int n = 1; n <= 2*t; n++)
    {
      termino = f(n);
      suma = suma + termino;
      cout << n << "\t" << termino << "\t" << suma <<endl;
    }
  
  return 0;
}

double f(int n)
{
  double funcion = 0.0;
  funcion = (2.0*n)/(2.*n+1);
  return funcion;
}
