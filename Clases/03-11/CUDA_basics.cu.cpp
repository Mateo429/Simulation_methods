#include <iostream>
#include<fstream>
#include <cmath>


#define Lx 16
#define Nx 8
const int Mx = (Lx+Nx-1)/Nx;

//----Programa del Device-------

//Kernel.

__global__ void AddTwoVectors(float *d_a, float *d_b, float *d_c){

   //¿Qué tarea me toca?

   int ix; ix = blockIdx.x*blockDim.x + threadIdx.x;
   dc[ix] = d_a[ix] + d_b[ix];

 }

//------Código del host-----

int main(){

  //Declarar todas las varaibles

  int ix;
  float h_a[Lx], h_b[Lx], h_c[Lx];

  //---en el device

  float *d_a; cudaMalloc((void**) &d_a, Lx*sizeof(float));
  float *d_b; cudaMalloc((void**) &d_b, Lx*sizeof(float));
  float *d_c; cudaMalloc((void**) &d_c, Lx*sizeof(float));

  //Inicializar variables

  for(ix =0; ix<Lx; ix++){
    h_a[ix] = ix; h_b[ix] = 2*ix;

  }
    
  //Enviar al device
  cudaMemcpy(d_a, h_a, Lx*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, h_b, Lx*sizeof(float), cudaMemcpyHostToDevice);

  //Correr en el device

  dim3 ThreadsPerBlock(Nx,0,0);
  dim3 BlocksPerGrid(Mx,0,0);
  AddTwoVectors<<<BlocksPerGrid, ThreadsPerBlock>>>(d_a, d_b, d_c);

  //Devolver el resultado al host
  
  cudaMemcpy(d_c, h_c, Lx*sizeof(float), cudaMemcpyDeviceToHost);

  //Imprimir resultados
  for(ix=0; ix<Lx; ix++)
    std::cout<<h_c[ix]<<endl;

  //Liberar la memoria
  cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
 


  return 0;
}


