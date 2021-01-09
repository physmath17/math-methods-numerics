/* Ising 3D code for MM2 course 2020 by Sayantan Sharma*/
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<complex.h>
#include<stdlib.h>
#include"ran3.h"
int nx,ny;
int main(int argc, char *argv[])
{
  int i,j,cutoff;
  int vol,sweep,maxsweep;
  double energy,de,r;
  double *sp,*spnew;
  double betaJ,M,Msq,Mav,H,sumH;
  double c_t,c0;
  long seed=-1.0;
  FILE *pt1;
  FILE *pt2;
  int iup(int ,int );
  int idn(int,int);
  //pt1=fopen(argv[4],"w");
  //pt2=fopen("results","r");
  //Lattice size
  nx=atoi(argv[1]);
  ny=atoi(argv[2]);
  //Value of inverse temperature*coupling
  betaJ=atof(argv[3]);
 
  //Total volume
   vol=nx*ny;
  
  //Total Monte Carlo sweeps
  maxsweep=300000;
 
  sp=malloc(vol*sizeof(double));
  spnew=malloc(vol*sizeof(double));
  
  
  
  sweep=0;
  H=0.0;
  
  j=0;
  Mav=0.0;
  M=0.0;
  sumH=0.0;
  //Starting configuration generation
  for(i=0;i<vol;i++){
    //hot start
    //if(ran3(&seed)>0.5) sp[i]=1.0 else sp[i]=-1.0;

    //cold start
    sp[i]=1.0;

    M+=sp[i];
  }
  
  while(sweep<maxsweep)
    {
      Msq=0.0;
      
      H=0.0;
      for(i=0;i<vol;i++)
	{
	  
      spnew[i]=-sp[i];
      M-=sp[i];
      //Compute the change in the energy due to spin flip
      //Note that 0 is the x-direction and 1 is along y-direction
      de=-(spnew[i]-sp[i])*(sp[iup(i,0)]+sp[idn(i,0)]+sp[iup(i,1)]+sp[idn(i,1)]);
      r=ran3(&seed);
      //Metropolis aceept reject criterion
      if(exp(-betaJ*de)>=r) sp[i]=spnew[i];
      M+=sp[i];
      Msq+=sp[i]*sp[i];
      //for total energy
      H+=sp[i]*(sp[iup(i,0)]+sp[idn(i,0)]+sp[iup(i,1)]+sp[idn(i,1)]);
    }
  
  
  
  cutoff=1e5;
  
  if((sweep)%200==0&&sweep>cutoff)
    {
      
      Mav+=fabs(M);
      sumH+=H;
      j+=1;
    }
  
  //fprintf(pt1,"%d\t%g\n",sweep,Msq/vol-pow(Mav/vol,2));
  //printf("%d\t%g\n",sweep,Msq/vol-pow(Mav/vol,2));
  sweep+=1;
    }	 	 
  printf("The value of magnetization= %g ,energy= %g at betaJ= %g \n",Mav/(j*vol),sumH/(j*vol),betaJ);
  //Read and construct specific heat
  /*sumH=0.0;
  sumHsq=0.0;
  while(eof(pt2))
    {
      fscanf(pt2, "%d\t%g\t%g\n",&t,&c_t,&H);
      sumH+=H;
      sumHsq+=H*H;
    }
    Cv=(sumHsq/vol-pow(sumH/vol,2))/vol;*/
    
  //fclose(pt1);
  free(sp);
  free(spnew);
  return(0);
}
  
int iup(int n1,int n2)
{
  int x,y;
  y=(int)n1/nx;
  x=n1-y*nx;
  if(n2==0) return(n1+(x+1)%(nx)-x);
  if(n2==1) return (n1+nx*((y+1)%ny-y));
}

int idn(int n1,int n2)
{
    int x,y;
    y=(int)n1/nx;
    x=n1-y*nx;
    if(n2==0) return(n1+(x-1+nx)%nx-x);
    if(n2==1) return (n1+nx*((y-1+ny)%ny-y));
    
}  
