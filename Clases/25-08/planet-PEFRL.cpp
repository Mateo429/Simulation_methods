#include <iostream>
#include<cmath>
#include "vector.h"

class Body{
private:
  vector3D r, V, F; 
  double m, R;
public:

void Initialize(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,  double R0);
void Compute_force(double GM);
void Start(double dt);  
void move_r(double dt, double theta);
void move_v(double dt, double theta);
  double Getx(void){return r.x();}; 
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  

  

};

void Body::Initialize(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,  double R0){

r.load(x0, y0, z0);
V.load(Vx0, Vy0, Vz0);
m = m0;
R = R0;
}

void Body::Compute_force(double GM){

  double aux = GM*m*std::pow(r.norm2(), -1.5);
  F = -(aux)*r;

  
}

void Body::Start(double dt){

  V -= -1*F*(0.5/m)*dt;
}

void Body::move_r(double dt, double theta){


  r+=V*(dt*theta);

}

void Body::move_v(double dt, double theta){

  V += F*(dt*theta/m);
}

int main() {
Body Planet;
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coef1=(1-2*Lambda)/2;
const double Coef2=1-2*(Chi+Lambda);
double t, dt = 0.01;
double GM = 1;
double R = 1;
double Omega = std::sqrt(GM/std::pow(R,3));
double T = 2*M_PI/Omega;
 Planet.Initialize(R, 0, 0, 0, R*Omega/2, 0, 0.453, 1);
 
  
  for(t=0; t<T; t+=dt){
    std::cout<<Planet.Getx()<<" "<< Planet.Gety()<< " "<< " "<< Planet.Getz()<<std::endl;
    Planet.move_r(dt, Zeta);
    Planet.Compute_force(GM);
    Planet.move_v(dt, Coef1);
    Planet.move_r(dt, Chi);
    Planet.Compute_force(GM);
    Planet.move_v(dt, Lambda);
    Planet.move_r(dt, Coef2);
    Planet.Compute_force(GM);
    Planet.move_v(dt, Lambda);
    Planet.move_r(dt, Chi);
    Planet.Compute_force(GM);
    Planet.move_v(dt, Coef1);
    Planet.move_r(dt, Zeta); 
    

    
  }
  }

