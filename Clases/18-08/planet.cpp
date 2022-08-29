#include <iostream>
#include<cmath>

class Body{
private:
double x, y, Vx, Vy, Fx, Fy, m, R;
public:

void Initialize(double x0, double y0, double Vx0, double Vy0,double m0,  double R0);
void Compute_force(double GM, double R);
void Move(double dt);
double Getx(void){return x;}; 

double Gety(void){return y;};

};



void Body::Initialize(double x0, double y0, double Vx0, double Vy0, double m0,  double R0){

x=x0;
y=y0;
Vx= Vx0;
Vy = Vy0;
m = m0;
R = R0;
}

void Body::Compute_force(double GM, double R){
Fx = -GM*m*x/(std::pow(R, 3));
Fy = -GM*m*y/(std::pow(R, 3)); 
  
}

void Body::Move(double dt){
x += Vx*dt;  y += Vy*dt;

Vx += Fx*dt/m; Vy += Fy*dt/m;

  
  
}

int main() {
Body Planet;
double t, dt = 0.01;
double GM = 1;
double R = 1;
double Omega = std::sqrt(GM/std::pow(R,3));
double T = 2*M_PI/Omega;
Planet.Initialize(R, 0, 0, R*Omega/2, 0.453, 1);
  
  for(t=0; t<T; t+=dt){
    std::cout<<Planet.Getx()<<" "<< Planet.Gety()<< " "<<std::endl;
    Planet.Compute_force(GM, R);
    Planet.Move(dt);

    
  }
  }
