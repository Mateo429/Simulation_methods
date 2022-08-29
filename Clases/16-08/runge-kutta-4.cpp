#include<iostream>
#include<cmath>

double f(double t, double x);
void Runge_one_step(double &t0, double &x0, double dt);

int main(){

  double t,x;
  double dt = 0.1;
  for(t=0, x=1; t<2; ){

    std::cout<<t<<" "<< x << " "<<std::exp(t) <<std::endl;
      Runge_one_step(t,x,dt);
  }

  
  return 0;
}

double f(double t, double x){

  return x;
}

void Runge_one_step(double &t0, double &x0, double dt){

  double dx1, dx2, dx3, dx4;
  dx1 = f(t0,x0)*dt;
  dx2 = dt*f(t0+ 0.5*dt, x0 + 0.5*dx1);
  dx3 = dt*f(t0 + 0.5*dt, x0 + 0.5*dx2);
  dx4 = dt*f(t0 + dt, x0 + dx3);
  
  x0+=(dx1 + 2*dx2 + 2*dx3 + dx4)/6;
  t0+=dt;

}
