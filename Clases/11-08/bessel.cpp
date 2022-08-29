#include<iostream>
#include<cmath>
double f(double x);
double Simpson(double a, double b, int n);
double Bessel(double alpha, double x);
double Bisection(double a, double b, double alpha);


int main(){
  double alpha, x;
  alpha = 0;
  for (x=0; x<=10; x+=0.1){

    std::cout<<x<<" "<< Bessel(alpha, x)<<std::endl;
  }

  std::cout<<Bisection(2, 4, alpha);

  return 0; 
}
double f(double t, double x, double alpha){

  return cos(alpha*t-x*sin(t));
}

double Simpson(double a, double b, int n, double x, double alpha){
  double t = 0;
  double h = 0;
  double sum = 0;
  n*=2; h=(b-a)/n;
  for(int ii = 0; ii<=n; ii++){
    t = a + ii*h;
    if(ii==0 || ii == n){
      sum += f(t, x, alpha);
    }
    else if(ii%2==0){
      sum += 2*f(t, x, alpha);
      }
    else {
      sum += 4*f(t, x, alpha); 
    }
    
  
  }
   return sum*(h/3);

}

double Bessel(double alpha, double x){

  return (1/(M_PI))*Simpson(0, M_PI, 100, x, alpha);

}

double Bisection(double a, double b, double alpha){
  double m, Bessela, Besselm;
  double Errmax= 1e-7;
  Bessela=Bessel(alpha, a);
  while(b-a>Errmax){
    m = (b+a)/2; Besselm = Bessel(alpha, m);
    if(Bessela*Besselm > 0){
      a = m;
      Bessela = Besselm;
      
    }
    else{ b = m;}

   }
  return (a+b)/2;
}
