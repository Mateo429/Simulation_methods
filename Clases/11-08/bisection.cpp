#include<iostream>
#include<cmath>

double f(double x);
double bisection(double a, double b);

///////////////////////////////////////////////////////
int main(){

  std::cout<<bisection(2,4)<<std::endl;





  return 0;
}

///////////////////////////////////////////////////////


double f(double x){

  return sin(x)/x;

}
//////////////////////////////////////////////////////////

double bisection(double a, double b){
  double x, m, fa, fm;
  double Errmax= 1e-7;
  fa=f(a);
  while(b-a>Errmax){
    m = (b+a)/2; fm = f(m);
    if(fa*fm > 0){
      a = m;
      fa = fm;
      
    }
    else{ b = m;}

   }
  return (a+b)/2;
}


/////////////////////////////////////////////////////////////  
  


