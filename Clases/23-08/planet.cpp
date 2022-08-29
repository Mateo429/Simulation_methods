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
void Move(double dt);
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

void Body::Move(double dt){

r += V*dt; V += (dt/m)*F; 

  
  
}

int main() {
Body Planet;
double t, dt = 0.001;
double GM = 1;
double R = 1;
double Omega = std::sqrt(GM/std::pow(R,3));
double T = 2*M_PI/Omega;
 Planet.Initialize(R, 0, 0, 0, R*Omega, 0, 0.453, 1);
  
  for(t=0; t<T; t+=dt){
    std::cout<<Planet.Getx()<<" "<< Planet.Gety()<< " "<< " "<< Planet.Getz()<<std::endl;
    Planet.Compute_force(GM);
    Planet.Move(dt);

    
  }
  }
