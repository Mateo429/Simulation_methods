#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=600;
const int Ly=200;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  double w[Q];      //Weights 
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i, int ix, int iy);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Start(double rho0,double Jx0,double Jy0);
  void Print(const char * NombreArchivo);
  double Ccelda(int ix, int iy);
  };  

double LatticeBoltzmann::Ccelda(int ix, int iy){
  int ix0 = 100;
  double n1=1;
  double n2=2*n1;
  double m = -1/(tan(0.52));
  double b = (Ly - m*ix0);
  int x0 = 50;
  int y0 = 100;
  int R = 100;
  
  /* if(ix < ((iy-b)/m) and iy < m*ix +b){ return (C/n1); }
     else{ return (C/n1)*(1-0.5*tanh(ix-ix0)); }*/
  if(pow(ix-x0,2) + pow(iy-y0,2) < R*R and ix < x0+R){return C/n1;}
  else{return 0;}

  }
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}  
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i, int ix, int iy){
   double C_r = Ccelda(ix, iy);
  
  if(i>0)
    return 3*w[i]*(C_r*C_r*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*(1-3*(C_r*C_r)*(1-W0));
}  
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i);
	f[n0]=feq(rho0,Jx0,Jy0,i, ix, iy);
      }
}  
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i, ix, iy);
      }
    }  
}
void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0,ix1,iy1,n01;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10; omega=2*M_PI/lambda*C;
  //an oscillating source in the middle
  /* ix=20; iy=20;
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i, ix, iy); 
    }*/

   ix =0;
  
  for (iy =0; iy<Ly; iy++){
    //omega=2*M_PI/lambda*Ccelda(ix, iy);
    rho0=100*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i, ix, iy); 
    }
    }


  /* for (iy =0; iy<Ly; iy++){
     for(ix=0; ix<Lx; ix++){
  
       rho0=30*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
       if(iy == -ix + 30){    
        for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        fnew[n0]=feq(rho0,Jx0,Jy0,i, ix, iy); 
	}
       }
     }
     }*/
   
  
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<Lx/3;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}
//--------------- Global Functions ------------

int main(void){
  LatticeBoltzmann Waves;
  int t,tmax=500;
  double rho0=0,Jx0=0,Jy0=0;

  //Start
  Waves.Start(rho0,Jx0,Jy0);
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision();
    Waves.ImposeFields(t);
    Waves.Advection();
  }
  //Show
  Waves.Print("Waves2D.dat");
 
  return 0;
} 
