#include<iostream>
#include<cmath>

double f(double t, double x);
void Euler_one_step(double &t, double &x, double dt);

int main(){

  double t,x;
  double dt = 0.01;
  for(t=0, x=1; t<2; ){

    std::cout<<t<<" "<< x << " "<<std::exp(t) <<std::endl;
      Euler_one_step(t,x,dt);
  }

  
  return 0;
}

double f(double t, double x){

  return x;
}

void Euler_one_step(double &t, double &x, double dt){

  double dx;
  dx=f(t,x)*dt;
  x+=dx;
  t+=dt;

}
