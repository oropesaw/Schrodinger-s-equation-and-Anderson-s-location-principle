///////////////////////////////////////////////////////////////////
//  William Carreras Oropesa. 17/9/2016                          //
//  Instituto Balseiro. FISCOM                                   //
//  Principio de localizacion de Anderson.                       //
///////////////////////////////////////////////////////////////////

# include <stdio.h>
# include <complex.h>
# include <math.h>
# include <time.h>
# include <assert.h>
# include <stdlib.h>

# define L 10.0
# define N  10000
# define V0 2e4
# define h  1
# define sigma 0.05
# define K0 (50.0*M_PI)
# define S0_2  0.5
# define X0 (L/2.0)
# define X00 (L/4.0)
# define ALPHA (I*h)
# define BETTA (I*1e-6*h)

typedef double real;
typedef double complex complejo;

typedef complejo(* funtion)(real);
typedef real (* potencial)(real);

typedef struct system
{
    funtion F;
    potencial V;
    real *P, dx,dt;
    complejo *Psi, *B, *C, *D;
}* clase;

// Funciones de estucturamiento y creacion del sistema.
clase Creacion (real dx0, real dt0, funtion F, potencial V);
void Free_Syst (clase pS);
void Initializ (clase pS);
void Iteration (clase pS);
void Prin_Syst(clase pS, FILE *gnuplotpipe);

// Funciones de calculo
complejo F (real x);
complejo P (real x);
real V_Free (real x);
real V_Randon (real x);
real Normalizate(clase pS);
real Norma_Psi (clase pS);
real M_P (clase pS);
real M_P_2 (clase pS);

int main(int argc, char const *argv[])
{
  real i;
  int j=0;
  FILE *feop = fopen("Anderson_Randon","w");
  FILE *pipe = popen("gnuplot", "w");
  assert(pipe);

  srand(time(NULL));

  clase pS=Creacion(0.001*h, 2e-6*h, P , V_Randon);
  Initializ (pS);

  for (i = 0; i < 0.1; i+=pS->dt)
  { 
     Prin_Syst(pS,pipe);
    fprintf(feop,"%lf\t%lf\t%lf\t%lf\n",i,cabs(pS->Psi[5000]*conj(pS->Psi[5000])),M_P(pS), sqrt(fabs(M_P_2(pS)-M_P(pS)*M_P(pS)))); 
    Iteration(pS);
    j++;
  }
  
  fclose(feop);
 fclose(pipe);
  Free_Syst(pS);
	return 0;
}

clase Creacion (real dx0, real dt0, funtion F, potencial V)
{
    clase pS = (clase)malloc(sizeof(struct system));
    assert(pS);
    
    pS->P = (real *)malloc(N*sizeof(real));
    assert(pS->P);

    pS->Psi = (complejo *)malloc(N*sizeof(complejo));
    assert(pS->Psi);

    pS->B = (complejo *)malloc(N*sizeof(complejo));
    assert(pS->B);
    
    pS->C = (complejo *)malloc(N*sizeof(complejo));
    assert(pS->C);

    pS->D = (complejo *)malloc(N*sizeof(complejo));
    assert(pS->D);

    pS->dt = dt0;
    pS->dx = dx0;
    pS->F = F;
    pS->V = V;

    return pS;
}

void Initializ (clase pS)
{
	int i,x;
    real norma;
    
    for (i = 0; i < N; ++i)
    pS->Psi[i] = pS->F(i*pS->dx);
    
    norma =Norma_Psi(pS);

    for (i = 0; i < N; ++i)
    pS->Psi[i] = norma*pS->Psi[i];

    pS->Psi[0] = pS->Psi[N-1]=0;
   
    for (i = 0; i < N; ++i)
    pS->P[i] = pS->V((i)*pS->dx);
    
   
    for(i= 0; i < N ; ++i)
    {
       pS->B[i] = 1+ALPHA+BETTA*pS->P[i];
       pS->C[i] = -0.5*ALPHA;
    }
}

void Iteration (clase pS)
{
	int i;
    
    pS->C[1]=-(0.5*ALPHA)/pS->B[1];
    for (i = 2; i < N-1; ++i)
    pS->C[i]=-(0.5*ALPHA)/(pS->B[i]+0.5*ALPHA*pS->C[i-1]);

    for ( i = 1; i < N-1; ++i)
    pS->D[i]=0.5*(ALPHA*pS->Psi[i-1])+(1.0-ALPHA-BETTA*pS->P[i])*pS->Psi[i]+0.5*(ALPHA*pS->Psi[i+1]);

    pS->D[1]=pS->D[1]/pS->B[1];
    for ( i = 2; i < N-1; ++i)
    pS->D[i]=(pS->D[i]+0.5*ALPHA*pS->D[i-1])/(pS->B[i]+0.5*ALPHA*pS->C[i-1]);
    
    pS->Psi[N-2]=pS->D[N-2];
    for ( i = N-3; i >0; i--)
    pS->Psi[i]=pS->D[i]-pS->C[i]*pS->Psi[i+1];
}

real Norma_Psi (clase pS)
{
  int i;
  real suma=0,norma;
  complejo temp;
  
  for (i= 1; i < N-1; ++i){
  temp=pS->Psi[i];
  suma+=cabs(temp*temp*pS->dx);
  }
  return 1.0/sqrt(suma);
}


real Normalizate(clase pS)
{
    real cont=0;
    complejo temp;
    int i;

    for (i = 1; i < N-1; ++i)
    {
       temp=pS->Psi[i];
       cont+=cabs(temp*temp*pS->dx);
    }
    return cont;
}

real V_Free (real x)
{
	return 0;
}

real V_Randon (real x)
{ 
   real m = ((1.0*rand())/RAND_MAX);
   return (V0/sqrt(0.001))*(2*m-1);
}

complejo F (real x)
{
  if(x==0.5*L) return (1.0/sqrt(0.001));
  else return 0;
}

real M_P (clase pS)
{
    int i;
   real Mean_position=0,aux;

   for (i = 1; i < N-1; i++)
   { 
        aux=cabs(pS->Psi[i]);
        Mean_position+=1.0*pS->dx*i*aux*aux*pS->dx;
   }
 
 return Mean_position;
}

real M_P_2 (clase pS)
{
   int i;
   real Mean_position=0,aux;

   for (i = 1; i < N-1; i++)
   { 
        aux=cabs(pS->Psi[i]);
        Mean_position+=1.0*pS->dx*i*pS->dx*i*aux*aux*pS->dx;
   }
 
 return Mean_position;
}

complejo P (real x)
{
   return cexp(ALPHA*K0*x)*exp(-((x-X00)*(x-X00))/(4*sigma*sigma));
}

void Prin_Syst(clase pS, FILE *gnuplotpipe)
{
  int i,x;

  fprintf(gnuplotpipe, "set title '{/=20 Evolucion temporal}'\n");
  fprintf(gnuplotpipe, "set tics \n");
  fprintf(gnuplotpipe, "set xrange[0:10]\n");
  fprintf(gnuplotpipe, "set ylabel  '{/=15 Densidad de Proba}'\n");
  fprintf(gnuplotpipe, "set xlabel  '{/=15 Posicion}'\n");
  fprintf(gnuplotpipe,  "plot '-' u 1:2  w l lw 2  t 'Norma=%lg'\n",Normalizate(pS));

  for (i=0;i<N;i++)
  {
   fprintf(gnuplotpipe, "%lf\t%lf\n", i*pS->dx,cabs(pS->Psi[i]*conj(pS->Psi[i])));
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

