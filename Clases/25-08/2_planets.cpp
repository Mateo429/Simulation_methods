#include <iostream>
#include<cmath>
#include "vector.h"

class Body;
class Colisionator;


class Body{
private:
  vector3D r, V, F; 
  double m, R;
public:

void Initialize(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,  double R0);
void erase_force(void);
void sum_force(vector3D F0);  
void Start(double dt);
void Draw(void);
void move_r(double dt, double theta);
void move_v(double dt, double theta);
  double Getx(void){return r.x();}; 
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  friend class Colisionator;
  

};

class Colisionator{

private:
  
public:

  void Compute_forces(Body *Planet, int N, double G);
  void Compute_force_between(Body & Planet1, Body & PLanet2, double G);
  

};

void Body::Initialize(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,  double R0){

r.load(x0, y0, z0);
V.load(Vx0, Vy0, Vz0);
m = m0;
R = R0;
}

void Body::erase_force(void){

F.load(0,0,0);
}

void Body::sum_force(vector3D F0){

  F += F0;

  
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


void Colisionator::Compute_forces(Body *Planet, int N, double G){
  int i, j;
  for (i=0; i<N; i++){
    Planet[i].erase_force();

  }

  for (i =0; i<N; i++)
    for (j = i + 1; j<N; j++)
      Compute_force_between(Planet[i], Planet[j], G);

}

void Colisionator::Compute_force_between(Body &Planet1, Body &Planet2, double G){

  vector3D r21, n, F1; double d21,  F;
  r21 = Planet2.r - Planet1.r; d21 = r21.norm(); n = r21/d21;
  F = G*Planet1.m*Planet2.m*std::pow(d21, -2.0);
  F1 =F*n; Planet1.sum_force(F1); Planet2.sum_force(-1*F1);

}

///////////////////////////////////////////////////////////////////////////////////////////
void Start_animation(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Balon.gif'"<<endl;

  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-20:20]"<<std::endl;
  std::cout<<"set yrange[-20:20]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;  
}

void Body::Draw(void){

  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";

}

void Start_square(void){

  std::cout<<"plot 0,0";
  
}

void End_square(void){

  std::cout<<std::endl;

}

int main() {
int N = 2;  
Body Planet[N];
Colisionator Newton;
 double G = 0.0001;
 double m0 = 10, m1 = 1, r = 11;
 double M = m0 + m1, x0 = -m1*r/M, x1 = m0*r/M;
 double omega = std::sqrt((G*M/(r*r*r))), T = 2*M_PI/omega, V0 = omega*x0, V1 = omega*x1;
 double t_draw , t_square = T/100;
 
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coef1=(1-2*Lambda)/2;
const double Coef2=1-2*(Chi+Zeta);
double t;
double dt = 0.001;

 Planet[0].Initialize(x0, 0, 0, 0, 0.5*V0, 0, m0, 1);
 Planet[1].Initialize(x1, 0, 0, 0, 0.5*V1, 0, m1, 0.5);

 Start_animation();
 
 
 for(t=0, t_draw=0; t<10*T; t+=dt, t_draw+=dt){
   // std::cout<<Planet[1].Getx()<<" "<< Planet[1].Gety()<< " "<< " "<< Planet[1].Getz()<<std::endl;
   if(t_draw>t_square){
     Start_square();
     for(int i =0; i<=N; i++) Planet[i].Draw();
     End_square();
     t_draw = 0;
     
   }


   for (int i = 0; i<N; i++) Planet[i].move_r(dt, Zeta);
    Newton.Compute_forces(Planet, N, G);
    for (int i = 0; i<N; i++) Planet[i].move_v(dt, Coef1);
    for (int i = 0; i<N; i++) Planet[i].move_r(dt, Chi);
    Newton.Compute_forces(Planet, N, G);
    for (int i = 0; i<N; i++) Planet[i].move_v(dt, Lambda);
    for (int i = 0; i<N; i++) Planet[i].move_r(dt, Coef2);
    Newton.Compute_forces(Planet, N, G);
    for (int i = 0; i<N; i++) Planet[i].move_v(dt, Lambda);
    for (int i = 0; i<N; i++) Planet[i].move_r(dt, Chi);
    Newton.Compute_forces(Planet, N, G);
    for (int i = 0; i<N; i++) Planet[i].move_v(dt, Coef1);
    for (int i = 0; i<N; i++) Planet[i].move_r(dt, Zeta); 
    

    
  }
  }

