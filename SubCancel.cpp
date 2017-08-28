#include <iostream>
#include <cmath>

using namespace std;

//t es el numero de terminos

int t = 30;
double s1r;
double s2r;
double s_1;
double s_2;
double s_3;

double S1(int N);
double S2a(int N);
double S2b(int N);
double S3(int N);

int main()
{
  std::cout.precision(16);  std::cout.setf(std::ios::scientific);

  for(int i = 1; i <= t; i++)
    {
      s_1 = S1(i);
      s_2 = S2b(i)-S2a(i);
      s_3 = S3(i);
	
      s1r = (abs(s_1-s_3))/(s_3);
      s2r = (abs(s_2-s_3))/(s_3);

      cout << i << "\t" << s1r << "\t" << "\t" << s2r <<endl;
      
    }
  //  cout << S1(t) <<endl;
  //  cout << S2(t) <<endl;
  //  cout << S3(t) <<endl;
  return 0;
}

double S1(int N)
{
  double suma = 0;
  for(int n = 1; n <= 2*N; n++)
    {
      suma = suma + pow(-1.0,n)*(n/(n+1.0));
    }
  
  return suma;
}

double S2a(int N)
{
  double suma1 = 0;
  
  for(int n = 1; n <= N; n++)
    {
      suma1 = suma1 + (2.0*n-1)/(2.*n);
    }

  return suma1;
}

double S2b(int N)
{
  double suma2 = 0;
  
  for(int k = 1; k <= N; k++)
    {
      suma2 = suma2 + (2.0*k)/(2.*k+1);
    }
  
  return suma2;
}

/*
double S2(int N)
{
  double suma;
  double suma1 = 0;
  double suma2 = 0;;
  
  for(int n = 1; n <= N; n++)
    {
      suma1 = suma1 + (2.0*n-1)/(2.*n);
    }
  
  for(int k = 1; k <= N; k++)
    {
      suma2 = suma1 + (2.0*k)/(2.*k+1);
    }

  suma = suma2 - suma1;
  
  return suma;
}
*/


double S3(int N)
{
  double suma = 0;
  for(int n = 1; n <= N; n++)
    {
      suma = suma + (1.0)/((2.*n)*(2.*n+1));
    }
  
  return suma;
}
