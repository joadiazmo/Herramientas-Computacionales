#include <iostream>
#include <cmath>

using namespace std;

int fact(int x)
{
  int resultado = 1;
  for(int i = x; i>=1; i--)
    {
      resultado=resultado*i;
    }
  return resultado;
}

int main()
{
  int numero;
  cin>>numero;
  cout << fact(numero) <<endl;
  return 0;
}
/*
double exp(double x, int N)
{
  for(int i = 0; i<=N; i++)
    {
      
    }
}
*/
