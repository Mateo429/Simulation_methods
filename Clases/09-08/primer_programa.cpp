 #include<iostream>
#include<cmath>
double f(double x);

int main(){
  const double ErrMax=1e-7;
  double x;
  double a = 2;
  double b = 4;
  double m;
  double fa= f(a);
  

  while (b-a>ErrMax){

    m = (a+b)/2;
    double fm = f(m);
    if (fa*fm>0){
      a=m;
      fa=fm;
    }
    else{b = m;}

  }
   std::cout<<"EL cero es "<<(a+b)/2;
  

  



  return 0;

}

double f(double x){

  return sin(x)/x;
}
