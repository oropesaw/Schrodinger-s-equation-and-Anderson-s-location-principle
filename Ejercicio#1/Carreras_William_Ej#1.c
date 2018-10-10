////////////////////////////////////////////////////////////////////////
//  William Carreras Oropesa. 17/9/2016                               //
//  Instituto Balseiro. FISCOM                                        //
//  Solucion de la ecuacion de Schroedinger para la particula libre   //
//  y potencial constante de barrera.                                 //
////////////////////////////////////////////////////////////////////////

# include <stdio.h>
# include <complex.h>
# include <math.h>
# include <assert.h>
# include <stdlib.h>

# define L 1.0
# define N 1000
# define X0 0.25
# define K0 50.0*M_PI
# define sigma 0.05
# define h 1
# define ALPHA I
# define BETTA (I*1e-6)
# define LIM 0.032
# define V0 (2500*M_PI*M_PI)

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

// Funciones de calculo del sistema.
complejo F (real x);
real V_Free (real x);
real V_Barrera(real x);
real Norma_Psi (clase pS);
real Normalizate (clase pS);           // Calcula la norma de la funcion de ondas.
real M_P (clase pS);                   // Valor medio del operador posicion.
real M_P_2 (clase pS);                 // Valor medio del cuadrado del operador posicion.
real M_M (clase pS);                   // Valor medio del operador momentun.
real M_M_2 (clase pS);                 // Valor medio del cuadrado del operador momentun.
real M_H (clase pS);                   // Valor medio del operador hamiltoniano.

int main(int argc, char const *argv[])
{
  real i;
  int j=0;
  FILE *fout = fopen("EJ_#1","w");
  FILE *pipe = popen("gnuplot", "w");
  assert(pipe);

  clase pS=Creacion(0.001*h, 2e-6*h, F, V_Free);
  Initializ (pS);
  
  for (i = 0; i < (50*L)/(2.0*h*K0 ); i+=pS->dt)
  {  
     
    if(j%10==0)  Prin_Syst(pS,pipe);
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i, Normalizate(pS),M_H(pS), M_P(pS), M_P_2(pS),
    sqrt(fabs(M_P_2(pS)-M_P(pS)*M_P(pS))), M_M(pS), M_M_2(pS), sqrt(fabs(M_M_2(pS)-M_M(pS)*M_M(pS))),h);
    Iteration(pS);
    j++;
  }

  fclose(pipe);
  fclose(fout);
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
	int i;
	real norma;

	for (i = 0; i < N; ++i)
    pS->Psi[i] = pS->F(i*pS->dx);
    
    
    norma = Norma_Psi(pS);

    for (i = 0; i < N; ++i)
    {
       pS->Psi[i] = norma*pS->Psi[i];
       pS->P[i] = pS->V(i*pS->dx);
    }
    
    pS->Psi[0] = pS->Psi[N-1]=0;
   
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
    for ( i = N-3; i >0; --i)
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

real M_H (clase pS)
{
  int i;
  real m_h=0;
  complejo g=ALPHA/(2e-6);

  for (i = 1; i < N; ++i)
  m_h+=creal((conj(pS->Psi[i]))*(h*I)*(g*pS->Psi[i-1]-2.0*g*pS->Psi[i]+g*pS->Psi[i+1]-(I/h)*pS->P[i]*pS->Psi[i])*pS->dx);
  return m_h;
}

real M_M (clase pS)
{
   int i;
   real M_m=0;
   
   for (i = 1; i < N-1; ++i)
   M_m+=creal(-conj(pS->Psi[i])*(I*h*pS->Psi[i+1]-I*h*pS->Psi[i]));
   
   return M_m;
}

real M_M_2 (clase pS)
{
   int i;
   real M_m=0;
   
   for (i = 1; i < N-1; ++i)
   M_m+=creal(-conj(pS->Psi[i])*(h*pS->Psi[i+1]-2.0*h*pS->Psi[i]+h*pS->Psi[i-1])/pS->dx);
   
   return M_m;
}

complejo F (real x)
{
   return cexp(ALPHA*K0*x)*exp(-((x-X0)*(x-X0))/(4*sigma*sigma));
}

real V_Free (real x)
{
	return 0;
}

real V_Barrera (real x)
{
	return ((fabs(x-2*X0)<LIM)? V0:0);
}



void Prin_Syst(clase pS, FILE *gnuplotpipe)
{
  int i,x;

  fprintf(gnuplotpipe, "set title '{/=20 Evolucion temporal}'\n");
  fprintf(gnuplotpipe, "set tics \n");
  fprintf(gnuplotpipe, "set yrange[0:4] ; set xrange[0:1]\n");
  fprintf(gnuplotpipe, "set ylabel  '{/=15 |Psi|^2}'\n");
  fprintf(gnuplotpipe, "set xlabel  '{/=15 Posicion}'\n");
  fprintf(gnuplotpipe,  "plot '-' u 1:2  w l lw 2  t 'Mean Position=%lg'\n", M_P(pS));

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
