///////////////////////////////////////////////////////////////////
//  William Carreras Oropesa. 17/9/2016                          //
//  Instituto Balseiro. FISCOM                                   //
//  Solucion de la ecuacion de Schroedinger para el doble pozo.  //
///////////////////////////////////////////////////////////////////


# include <stdio.h>
# include <complex.h>
# include <math.h>
# include <assert.h>
# include <stdlib.h>

# define L 20
# define N 2000
# define S0_2  0.5
# define X0 -1.0


typedef double real;
typedef double complex complejo;

typedef complejo(* funtion)(real);
typedef real (* potencial)(real);

typedef struct system
{
    funtion F;
    potencial V;
    real *P, dx,dt,h;
    complejo *Psi, *B, *C, *D;
}* clase;


// Funciones de estucturamiento y creacion del sistema.
clase Creacion (real dx0, real dt0, funtion F, potencial V, real h);
void Free_Syst (clase pS);
void Initializ (clase pS);
void Iteration (clase pS);
void Prin_Syst(clase pS, FILE *gnuplotpipe);

// Funciones de calculo del sistema.
complejo F (real x);
real V (real x);
real Norma_Psi (clase pS);
real Normalizate(clase pS);
real M_P (clase pS);
real M_P_2 (clase pS); 
real Probabilidad (clase pS, real x0, real x);

int main(int argc, char const *argv[])
{
  real i, h ;
  int j=0;
  FILE *fout = fopen("Archivo","w");
  FILE *pipe = popen("gnuplot", "w");
  assert(pipe);
  
  for(h = 0.2; h <= 3; h+=0.2){
  clase pS=Creacion(0.01*h, 4e-4*h, F, V, h);
  Initializ (pS);
  
  for (i = 0; i < 10; i+=pS->dt)
  {  
     
    if(j%40==0)  Prin_Syst(pS,pipe);
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i, Normalizate(pS), M_P(pS), sqrt(fabs(M_P_2(pS)-M_P(pS)*M_P(pS))), Probabilidad(pS, 0, L), h);
    Iteration(pS);
    j++;
  }
  Free_Syst(pS);
 }
 
 pclose(pipe);
 fclose(fout);
	return 0;
}

clase Creacion (real dx0, real dt0, funtion F, potencial V, real h)
{
    clase pS = (clase)malloc(sizeof(struct system));
    assert(pS);
    
    pS->P = (real *)malloc(2*N*sizeof(real));
    assert(pS->P);

    pS->Psi = (complejo *)malloc(2*N*sizeof(complejo));
    assert(pS->Psi);

    pS->B = (complejo *)malloc(2*N*sizeof(complejo));
    assert(pS->B);
    
    pS->C = (complejo *)malloc(2*N*sizeof(complejo));
    assert(pS->C);

    pS->D = (complejo *)malloc(2*N*sizeof(complejo));
    assert(pS->D);

    pS->dt = dt0;
    pS->dx = dx0;
    pS->h = h;
    pS->F = F;
    pS->V = V;

    return pS;
}

void Initializ (clase pS)
{
	int i;
	real norma;
  
  complejo ALPHA = 2*I*pS->h;
  complejo BETTA = 0.5*I*4e-4/pS->h;

	for (i = 0; i < 2*N; ++i)
    pS->Psi[i] = pS->F((i-N)*pS->dx);
    
    
    norma = Norma_Psi(pS);

    for (i = 0; i < 2*N; ++i)
    {
       pS->Psi[i] = norma*pS->Psi[i];
       pS->P[i] = pS->V((i-N)*pS->dx);
    }
    
    pS->Psi[0] = pS->Psi[2*N-1]=0;
   
    for(i= 0; i < 2*N ; ++i)
    {
       pS->B[i] = 1+ALPHA+BETTA*pS->P[i];
       pS->C[i] = -0.5*ALPHA;
    }
}

void Iteration (clase pS)
{
	int i;
    
    complejo ALPHA = 2*I*pS->h;
    complejo BETTA = 0.5*I*4e-4/pS->h;

    pS->C[1]=-(0.5*ALPHA)/pS->B[1];
    for (i = 2; i < 2*N-1; ++i)
    pS->C[i]=-(0.5*ALPHA)/(pS->B[i]+0.5*ALPHA*pS->C[i-1]);

    for ( i = 1; i < 2*N-1; ++i)
    pS->D[i]=0.5*(ALPHA*pS->Psi[i-1])+(1.0-ALPHA-BETTA*pS->P[i])*pS->Psi[i]+0.5*(ALPHA*pS->Psi[i+1]);

    pS->D[1]=pS->D[1]/pS->B[1];
    for ( i = 2; i < 2*N-1; ++i)
    pS->D[i]=(pS->D[i]+0.5*ALPHA*pS->D[i-1])/(pS->B[i]+0.5*ALPHA*pS->C[i-1]);
    
    pS->Psi[2*N-2]=pS->D[2*N-2];
    for ( i = 2*N-3; i >0; --i)
    pS->Psi[i]=pS->D[i]-pS->C[i]*pS->Psi[i+1];
}

real Norma_Psi (clase pS)
{
  int i;
  real suma=0,norma;
  complejo temp;
  
  for (i= 1; i < 2*N-1; ++i){
  temp=pS->Psi[i];
  suma+=cabs(temp*temp*pS->dx);
  }
  return 1.0/sqrt(suma);
}

complejo F (real x)
{
	return exp(-((x-X0)*(x-X0))/(4*S0_2));
}

real V (real x)
{   
	real a=x*x;
	return -0.25*a+0.125*a*a;
}

real Normalizate(clase pS)
{
    real cont=0;
    complejo temp;
    int i;

    for (i = 1; i < 2*N-1; ++i)
    {
       temp=pS->Psi[i];
       cont+=cabs(temp*temp*pS->dx);
    }
    return cont;
}

real M_P(clase pS)
{
    int i;
   real Mean_position=0,aux;

   for (i = 1; i < 2*N-1; i++)
   { 
        aux=cabs(pS->Psi[i]);
        Mean_position+=1.0*pS->dx*(i-N)*aux*aux*pS->dx;
   }
 
 return Mean_position;
}

real M_P_2 (clase pS)
{
   int i;
   real Mean_position=0,aux;

   for (i = 1; i < 2*N-1; i++)
   { 
        aux=cabs(pS->Psi[i]);
        Mean_position+=1.0*pS->dx*(i-N)*pS->dx*(i-N)*aux*aux*pS->dx;
   }
 
 return Mean_position;
}

real Probabilidad (clase pS, real x0, real x)
{
	int i;
	real prob=0;
	complejo temp;
    
	for (i =floor(x0/pS->dx) ; i < floor(x/pS->dx) ; ++i)
	{ 
		temp=pS->Psi[i+N];
    prob+=creal(temp*conj(temp)*pS->dx);
	}
    
  return prob;
}

void Prin_Syst(clase pS, FILE *gnuplotpipe)
{
  int i,x;

  fprintf(gnuplotpipe, "set title '{/=20 Evolucion temporal}'\n");
  fprintf(gnuplotpipe, "set tics \n");
  fprintf(gnuplotpipe, "set yrange[0:1] ; set xrange[-20:20]\n");
  fprintf(gnuplotpipe, "set ylabel  '{/=15 Densidad de Proba}'\n");
  fprintf(gnuplotpipe, "set xlabel  '{/=15 Posicion}'\n");
  fprintf(gnuplotpipe,  "plot '-' u 1:2  w l lw 2  t 'Norma=%lg'\n", Normalizate(pS));

  for (i=0;i<2*N;i++)
  {
   fprintf(gnuplotpipe, "%lf\t%lf\n", (i-N)*pS->dx,cabs(pS->Psi[i]*conj(pS->Psi[i])));
  }

  fprintf(gnuplotpipe, "e\n" );
	fflush(gnuplotpipe);
}

void Free_Syst (clase pS)
{
  free(pS->B);
  free(pS->C);
  free(pS->D);
  free(pS->Psi);
  free(pS->P);
  free(pS);
}