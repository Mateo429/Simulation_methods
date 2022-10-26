#include<iostream>
#include <cmath>
#include <fstream>
using namespace std;

const int Lx = 128;
const int Ly = 128;
const int Q = 5;
const double W0 = 1.0/3;
const double C = 0.5;
const double C2 = C*C;
const double AUX0 = 1-3*C2*(1-W0);
const double tau = 0.5;
const double Utau = 1.0/tau;
const double Umtau = 1-Utau;
//--------------------- Clase LatticeGas ------------


class LatticeBoltzmann{
private:
  double w[Q]; //Pesos de las velocidades
  int Vx[Q], Vy[Q]; // Vectores de velocidad, separados por componentes
  double *f, *fnew;  // Funciones de distribución.

public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q + i;}; //Función que nos permite recorrer el Array de distribución.
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double feq(double rho0, double Jx0, double Jy0, int i);
  void Inicie(double rho0, double Jx0, double Jy0);
  void Colision(void);
  void ImponerCampos(int t);
  void Adveccion(void);
  void Imprimir(const char * Namefile);


 };

void LatticeBoltzmann::Imprimir(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix, iy;
  for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      rho0 = rho(ix, iy, true);
      MyFile << ix << " " << iy << " " << rho0 << endl; 

    }
    MyFile<<endl;
  }
  MyFile.close();
  }

void LatticeBoltzmann::Adveccion(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}

void LatticeBoltzmann::ImponerCampos(int t){

  int i, ix , iy, n0;
  double lambda, omega, rho0, Jx0, Jy0; lambda = 10; omega = 2*M_PI/lambda*C;
  //Se pondrá una fuente oscilante en el medio
  rho0 = 10*sin(omega*t); Jx0 = Jx(ix, iy, false); Jy0 = Jy(ix, iy, false); //Solo imponemos el campo de densidad.
  for(i =0; i<Q; i++){
    n0 = n(ix, iy, i);
    fnew[n0] = feq(rho0, Jx0, Jy0, i);
  }
}

void LatticeBoltzmann::Colision(void){

  int ix, iy, i, n0; double rho0, Jx0, Jy0;
  for (ix=0; ix<Lx; ix++)
    for (iy=0; iy<Ly; iy++){

      rho0= rho(ix, iy, false); Jx0=Jx(ix, iy, false); Jy0=Jy(ix, iy, false);
      for(i =0; i<Q; i++){
	n0 = n(ix, iy, i);
	fnew[n0] = Umtau*f[n0] + Utau*feq(rho0, Jx0, Jy0, i);
      }
    }

}
void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){

  int ix , iy, i, n0;
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      for (i=0; i<Q; i++){
        n0 = n(ix, iy, i);
	  f[n0] = feq(rho0, Jx0, Jy0, i);

      }

}

double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i){
  
  if(i>0) return 3*C2*rho0*w[i] + 3*w[i]*(Vx[i]*Jx0 + Vy[i]*Jy0);
  else
    return rho0*AUX0;

}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){

  double sum ; int i , n0;
  for (sum =0, i=0; i<Q; i++){
    n0 = n(ix, iy, i);
    if (UseNew) sum += fnew[n0]; else sum += f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
 double sum ; int i , n0;
  for (sum =0, i=0; i<Q; i++){
    n0 = n(ix, iy, i);
    if (UseNew) sum += Vx[i]*fnew[n0]; else sum += Vx[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){

 double sum ; int i , n0;
  for (sum =0, i=0; i<Q; i++){
    n0 = n(ix, iy, i);
      if (UseNew) sum += Vy[i]*fnew[n0]; else sum += Vy[i]*f[n0];
  }
  return sum;
}


 LatticeBoltzmann::LatticeBoltzmann(){ //Método constructor de la clase, se agregan los pesos, las velocidades, y se crean las "matrices" de dist.
  //Se inicilizan los pesos
  w[0]= W0;  w[1]=w[2]=w[3]=w[4]=(1-W0)/4;
  //Se inicializan las velocidades
  Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
  Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
  //Se crean los arreglos dinámicos
  int ArraySize = Lx*Ly*Q;
  f = new double [ArraySize]; fnew = new double [ArraySize]; //Se pide memoria al heap
}

LatticeBoltzmann::~LatticeBoltzmann(){

  delete [] f; delete [] fnew; //Se libera la memoria.
}



//---------------------  Programa Principal ------------

int main(void){
  LatticeBoltzmann Ondas;
  double t, tmax= 1000;
  double rho0, Jx0, Jy0;
  Ondas.Inicie(rho0, Jx0, Jy0);
  for(t=0; t<tmax; t++){
    Ondas.Colision();
    Ondas.ImponerCampos(t);
    Ondas.Adveccion();

  }
    
  // Ondas.Imprimir("Ondas.dat");
  
  return 0;
}
