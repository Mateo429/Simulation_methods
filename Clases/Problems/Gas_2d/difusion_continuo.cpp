#include<iostream>
#include <cmath>
using namespace std;

const int Lx=1024;
const double p=0.5;

const int Q=2;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 (derecha) i=1 (izquierda)
  double  f[Lx][Q], fnew[Lx][Q]; // n[ix][i]
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma);
  void Show(void);
  void GrafiqueRho(void);
  void Colisione();
  void Adveccione(void);
  double rho(int ix); //Inline
  double Varianza(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;  V[1]=-1;
}


void LatticeGas::Inicie(int N,double mu,double sigma){
  for (int ix=0; ix<Lx; ix++){
    double rho0 = (N/sigma*sqrt(2*M_PI))*exp(-0.5*pow(ix-mu, 2));
    for(int i =0; i<Q; i++)
      f[ix][i]= rho0/Q;
   
  }
}

double LatticeGas::rho(int ix){

  double suma; int i;

  for(suma=0, i=0; i<Q; i++)
    suma+=f[ix][i];
  return suma;
}

void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    std::cout<<ix<<" "<<rho(ix)<<std::endl;
}
void LatticeGas::Colisione(){
  int ix, i, j;
  for(ix=0; ix<Lx; ix++)
    for(i=0; i<Q; i++){
      j = (i+1)%Q;
	fnew[ix][i] = f[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);   
    }
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int i=0;i<Q;i++)
      f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i];
}

double LatticeGas::Varianza(void){
  int ix; double N,Xprom,Sigma2;
  //Calcular N
  for(N=0,ix=0;ix<Lx;ix++)
    N+=rho(ix);
  //Calcular Xprom
  for(Xprom=0,ix=0;ix<Lx;ix++)
    Xprom+=ix*rho(ix);
  Xprom/=N;
  //Calcular Sigma2
  for(Sigma2=0,ix=0;ix<Lx;ix++)
    Sigma2+=pow(ix-Xprom,2.0)*rho(ix);
  Sigma2/=(N-1);
  
  return Sigma2;
}


//---------------------  Programa Principal ------------

int main(void){
  LatticeGas Difusion;
  int N=400; double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;

  Difusion.Inicie(N,mu,sigma);
  for(t=0;t<tmax;t++){
    // cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }
  Difusion.GrafiqueRho();
  
  return 0;
}
