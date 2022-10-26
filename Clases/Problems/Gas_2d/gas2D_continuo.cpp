#include<iostream>
#include<cmath>
using namespace std;


const int Lx = 300;
const int Ly = 300;
const double p = 0.25;
const double p0 = 0.25;
const int N = 2400;
const int Q = 4;

class LatticeGas2D{
private:
  int Vx[Q];
  int Vy[Q];
  double f[Lx*Ly*Q], fnew[Lx*Ly*Q];


public:
  LatticeGas2D(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
  void Inicie(int N, double mux, double muy, double sigma);
  double rho(int ix, int iy);
  void Colisione();
  void Adveccione(void);
  double Varianza(void);
  void GrafiqueRho(void);
  int getf(int ix, int iy, int i){return f[n(ix, iy, i)];};
  void show(void);
  
};

void LatticeGas2D::show(void){
  int k;
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++)
      cout<<getf(ix, 3 , i)<<" ";
    cout<<endl;
  }
  cout<<endl;

}
LatticeGas2D::LatticeGas2D(void){

    Vx[0]=1;  Vx[1]=0;  Vx[2]=-1; Vx[3]=0;
    Vy[0]=0;  Vy[1]=1;  Vy[2]=0;  Vy[3]=-1;

}

void LatticeGas2D::Inicie(int N, double mux, double muy, double sigma){
int ix, iy, i, k;
 for(ix=0; ix<Lx; ix++)
   for(iy=0; iy<Ly; iy++){
    double rho0x = N*(1/sigma*sqrt(2*M_PI))*exp(-0.5*pow(ix-mux, 2));
    double rho0y = N*(1/sigma*sqrt(2*M_PI))*exp(-0.5*pow(iy-muy, 2));
    for(int i =0; i<Q; i++)
      f[n(ix, iy, i)]=(rho0x)*(rho0y)/Q;
    

   }
 }




double LatticeGas2D::rho(int ix, int iy){

  double sum = 0;
  int k;
  for(int i=0; i<Q; i++){
    k = n(ix, iy, i);
    sum += f[k];
    
  }
  return sum;

}


void LatticeGas2D::Colisione(){
  int j, k, l;

  for(int ix=0;ix<Lx;ix++) 
     for(int iy=0;iy<Ly;iy++)
       for(int i=0;i<Q;i++) {
	 j=(i+1+Q)%Q;
	 k=(i-1+Q)%Q;
	 l=(i+2+Q)%Q;
         fnew[n(ix,iy,i)] =p0*f[n(ix,iy,i)]+p*(f[n(ix,iy,j)]+f[n(ix,iy,k)])+(1-2*p-p0)*f[n(ix,iy,l)];       	 
   }
}

void LatticeGas2D::Adveccione(void){
int ix,iy,i,ixnext,iynext,n0,n0next;

 for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ 
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0];
      }
}


double LatticeGas2D::Varianza(void){
int ix, iy; double N, Xprom, Yprom, Sigma2x, Sigma2y, Sigma2;
 double rhox, rhoy;
  for(N=0,ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      N+=rho(ix, iy);

  Xprom=0;
  Yprom=0;
  for(ix = 0; ix<Lx; ix++){
     
        for(iy = 0; iy<Ly; iy++)
            rhox += ix*rho(ix, iy);
        Xprom +=rhox;
        rhox = 0;
    }
  Xprom/=N;

  for(iy = 0; iy<Ly; iy++){
     
        for(ix = 0; ix<Lx; ix++)
            rhoy += iy*rho(ix, iy);
        Yprom +=rhoy;
        rhoy = 0;
    }
  Yprom/=N;

    Sigma2x = 0;
    Sigma2y = 0;

    for(ix = 0; ix<Lx; ix++){
     
        for(iy = 0; iy<Ly; iy++)
            rhox += rho(ix, iy);
        Sigma2x += pow((Xprom - ix),2)*rhox;
        rhox = 0;
    }
    Sigma2x/=(N-1);

    for(iy = 0; iy<Ly; iy++){
        
        for(ix = 0; ix<Lx; ix++)
            rhoy += rho(ix, iy);
        Sigma2y += pow((Yprom - iy),2)*rhoy;
        rhoy = 0;
    }
    Sigma2y/=(N-1);

    return Sigma2x+Sigma2y;
 }

void LatticeGas2D::GrafiqueRho(void){

for(int ix=0;ix<Lx;ix++)
   for(int iy =0; iy<Ly;iy++)

     std::cout<< ix << " " << iy << " " << rho(ix, iy)/(N*N) <<endl; 

}

//////////////////////////////////////////////////////////////

int main(){

  LatticeGas2D Difusion;
  double mux = Lx/2, muy=Ly/2, sigma = 16;
  int t, tmax=300;
  Difusion.Inicie(N, mux, muy, sigma);
  // Difusion.GrafiqueRho();
  // Difusion.show();
  for(t=0; t<tmax; t++){
    //u = ran64x.r();
    //cout<<t<<" "<< Difusion.getn(Lx/2,Ly/2, 1) <<endl;
    // cout<<t<<" "<< Difusion.rho(Lx/2+2,Ly/2) <<endl;
    //cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
    // Difusion.show();
  }

  Difusion.GrafiqueRho();



  return 0;
}

