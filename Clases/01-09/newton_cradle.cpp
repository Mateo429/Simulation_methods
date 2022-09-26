#include <iostream>
#include <cmath>


using namespace std;//Constantes globales

const int N=3;
const double g=980;
const double k = 1e7;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  // vector3D r,V,F;  double m,R;
  double theta, omega;
  double tau;
  double m, R, l, x0, I;
public:
  void Inicie(double theta0, double omega0, double m0, double R0, double l0, double x00);
  void Borretorque(void){tau=0;};
  void Sumetorque(double tau0){tau+=tau0;};
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return x0 + l*std::sin(theta);}; //Inline
  double Gety(void){return -l*std::cos(theta);};
  double Gettau(void){return tau;};
  friend class Colisionador;
};


void Cuerpo::Inicie(double theta0,double omega0, double m0,double R0, double l0, double x00){
  //r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
  theta = theta0; omega = omega0; m = m0; l = l0; x0 = x00; I=m*l*l; R = R0;
}
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(dt*coeficiente);
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<theta<<"), -"<<l/7<<"*t*cos("<<theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeTorques(Cuerpo * Planeta);
  void CalculeTorqueEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);    
};

void Colisionador::CalculeTorques(Cuerpo * Planeta){
  int i,j; double tau_0;
  //Borrar fuerzas
  for(i=0;i<N;i++){
  Planeta[i].Borretorque();
  tau_0 = -Planeta[i].l*Planeta[i].m*g*std::sin(Planeta[i].theta);
  Planeta[i].Sumetorque(tau_0);}
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=N-1;i>0;i--)
      CalculeTorqueEntre(Planeta[i],Planeta[i-1]);
}
void Colisionador::CalculeTorqueEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  double s = (Planeta2.Getx() + Planeta2.R)-(Planeta1.Getx() -  Planeta1.R); double F =0; 
  if (s>0) F = k*pow(s, 1.5);
  Planeta1.Sumetorque(F*Planeta1.l); Planeta2.Sumetorque(-F*Planeta2.l);

}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'DosPlanetas.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:25]"<<endl;
  cout<<"set yrange[-15:2]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m0=50, l0 = 12, r0=1.5, theta0 = -0.2, omega0=10;
  double T=2*M_PI*std::sqrt(l0/g);
  double t,tmax=5*T,dt=1e-5;
  double tdibujo,tcuadro=T/500;
  int i;
  
  //---------------(x0,y0,z0,Vx0,Vy0,Vz,m0,R0)
  Planeta[0].Inicie(-0.40, 0, m0, r0, l0, 0);

  for(i =1; i<N; i++){

    Planeta[i].Inicie(0,0, m0, r0, l0, 2*r0*i);
  }
  
  
  InicieAnimacion();
  
  for(t=0,tdibujo=0; t<10*tmax; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) Planeta[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }         
    
    // cout<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<endl;
    // Mover por PEFRL
    for(i=0;i<N;i++) Planeta[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTorques(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_theta(dt,Coeficiente2);
    Newton.CalculeTorques(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_theta(dt,Zeta);   
  }
  
  return 0;
}
