#include<iostream>
#include<cmath>

double f1(double t, double x1, double x2);
double f2(double t, double x1, double x2);
void Runge_one_step(double &t0, double &x10, double &x20, double dt);

int main(){

  double t,x1, x2;
  double dt = 0.1;
  for(t=0, x1=1, x2=0; t<7; ){

    std::cout<<t<<" "<< x1 << " " << x2 <<std::endl;
    Runge_one_step(t,x1,x2,dt);
  }

  
  return 0;
}

double f1(double t, double x1, double x2){

  double omega = 3;
  
  return -omega*omega*x2;
}

double f2(double t, double x1, double x2){

  return x1;
}

void Runge_one_step(double &t0, double &x10, double &x20, double dt){

  double dx11, dx21, dx31, dx41;
  double dx12, dx22, dx32, dx42;


  dx11 = f1(t0,x10,x20)*dt;                                  dx12 = f2(t0,x10,x20)*dt;

  dx21 = dt*f1(t0+ 0.5*dt, x10 + 0.5*dx11, x20 + 0.5*dx12);  dx22 = dt*f2(t0+ 0.5*dt, x10 + 0.5*dx11, x20 + 0.5*dx12);

  dx31 = dt*f1(t0+ 0.5*dt, x10 + 0.5*dx21, x20 + 0.5*dx22);  dx32 = dt*f2(t0+ 0.5*dt, x10 + 0.5*dx21, x20 + 0.5*dx22);

  dx41 = dt*f1(t0+ 0.5*dt, x10 + dx31, x20 + dx32);          dx42 = dt*f2(t0+ 0.5*dt, x10 + dx31, x20 + dx32);  

  
  x10+=(dx11 + 2*dx21 + 2*dx31 + dx41)/6;
  x20+=(dx12 + 2*dx22 + 2*dx32 + dx42)/6;
  t0+=dt;

}
