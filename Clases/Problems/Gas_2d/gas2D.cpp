#include<iostream>
#include "Random64.h"

const int Lx = 256;
const int Ly = 256;
const double p = 0.25;
const double p0 = 0.25;

const int Q = 4;

class LatticeGas2D{
private:
  int Vx[Q];
  int Vy[Q];
  int n[Lx*Ly*Q], nnew[Lx*Ly*Q];


public:
  LatticeGas2D(void);
  int indice(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
  void Inicie(int N, double mux, double muy, double sigma, Crandom & ran64);
  void Borrese(void);
  int rho(int ix, int iy);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double Varianza(void);
  void GrafiqueRho(void);
  int getn(int ix, int iy, int i){return n[indice(ix, iy, i)];};
  void show(void);
  
};

void LatticeGas2D::show(void){
  int k;
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++)
      cout<<getn(ix, 3 , i)<<" ";
    cout<<endl;
  }
  cout<<endl;

}
LatticeGas2D::LatticeGas2D(void){

    Vx[0]=1;  Vx[1]=0;  Vx[2]=-1; Vx[3]=0;
    Vy[0]=0;  Vy[1]=1;  Vy[2]=0;  Vy[3]=-1;

}

void LatticeGas2D::Inicie(int N, double mux, double muy, double sigma, Crandom &ran64){

  int ix, iy, i, k;
  while (N>0){
    ix = (int) ran64.gauss(mux, sigma);
    iy = (int) ran64.gauss(muy, sigma);
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1;
    if(iy<0) iy=0; if(iy>Ly-1) ix=Ly-1;
    i=(int) Q*ran64.r();
    k = indice(ix, iy, i);
    if(n[k]==0){n[k]=1; N--;}
   }
 }

void LatticeGas2D::Borrese(void){
  for(int ix=0;ix<Lx;ix++)
   for(int iy=0; iy<Ly; iy++)
     for(int i=0;i<Q;i++)
       n[indice(ix, iy, i)]=0;
}


int LatticeGas2D::rho(int ix, int iy){

  int sum = 0;
  int k;
  for(int i=0; i<Q; i++){
    k = indice(ix, iy, i);
    sum += n[k];
    
  }
  return sum;

}


void LatticeGas2D::Colisione(Crandom & ran64){
  int k; double u;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      u = ran64.r();
      k=indice(ix,iy,0);
    
     if(p0<u<p0+p){

     nnew[k] = n[k+3]; nnew[k+1] = n[k]; nnew[k+2] = n[k+1]; nnew[k+3] = n[k+2];
     }

     else if(p0+p<u<p0+2*p){

     nnew[k] = n[k+1]; nnew[k+1] = n[k+2]; nnew[k+2] = n[k+3]; nnew[k+3] = n[k];
     }

     else if(p0+2*p<u<1){

     nnew[k] = n[k+2]; nnew[k+1] = n[k+3]; nnew[k+2] = n[k]; nnew[k+3] = n[k+1];
     }

     else{

     nnew[k] = n[k]; nnew[k+1] = n[k+1]; nnew[k+2] = n[k+2]; nnew[k+3] = n[k+3];

     }
  
 }
}
}

void LatticeGas2D::Adveccione(void){
int ix,iy,i,ixnext,iynext,n0,n0next;

 for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ 
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=indice(ix,iy,i); n0next=indice(ixnext,iynext,i);
	n[n0next]=nnew[n0];
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
  /* for(ix=0;ix<Lx;ix++)
        for(iy=0; iy<Ly; iy++){
          Xprom+=ix*rho(ix, iy);
          Yprom+=iy*rho(ix, iy);}
    Xprom/=N;
    Yprom/=N;*/
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
    Sigma2x/=N;

    for(iy = 0; iy<Ly; iy++){
        
        for(ix = 0; ix<Lx; ix++)
            rhoy += rho(ix, iy);
        Sigma2y += pow((Yprom - iy),2)*rhoy;
        rhoy = 0;
    }
    Sigma2y/=N;

    return Sigma2x+Sigma2y;
   
  
  

}

void LatticeGas2D::GrafiqueRho(void){

for(int ix=0;ix<Lx;ix++)
   for(int iy =0; iy<Ly;iy++)

     std::cout<< ix << " " << iy << " " << rho(ix, iy) <<endl; 

}

//////////////////////////////////////////////////////////////

int main(){

  LatticeGas2D Difusion;
  Crandom ran64(0);
  int N = 2400; double mux = Lx/2, muy=Ly/2, sigma = 16;
  int t, tmax=300;
  double u;
  Difusion.Borrese();
  Difusion.Inicie(N, mux, muy, sigma, ran64);
  // Difusion.GrafiqueRho();
  // Difusion.show();
  for(t=0; t<tmax; t++){
    //u = ran64x.r();
    //cout<<t<<" "<< Difusion.getn(Lx/2,Ly/2, 1) <<endl;
    // cout<<t<<" "<< Difusion.rho(Lx/2+2,Ly/2) <<endl;
    cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    // Difusion.show();
  }

  //Difusion.GrafiqueRho();



  return 0;
}

