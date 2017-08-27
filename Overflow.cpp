#include <iostream>
#include <cstdlib>

using namespace std;

float under = 1.0;
float over = 1.0;
int N = 200;

int main()
{
  for(int i =1; i < N; i++)
    {
      under = under/2;
      over = over*2;
      
      cout << i << "\t" << under << "\t"  << over <<endl;
    } 
  return 0;
}
