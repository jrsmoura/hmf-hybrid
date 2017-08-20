// HMFNonIdentical.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ran2.h"
#include "Init.h"
#include "Energy.h"
#include"Force.h"
#include "Integration.h"
#include "Virial.h"

double Virial(long, double *, double *);
double TimeAverage(long, double *, double *);

using namespace std;

#define pi	3.14159265359
#define dpi	6.28318530718

int main(int argc, char* argv[])
{
  long n, seed, idum, Np;
  double p0, r0, alpha;
  double energKin, energPot, magX, magY, energ, energ0, error;
  double time, finalTime, timeStep, timeOutput, timeCount;
  double pm, dpSize, tMean, npMean, var, pMean, Nt, pp, statMoment4;

  FILE *init = fopen("./initialPhase.dat", "w");
  FILE *enrg = fopen("./energy.dat", "w");
  FILE *finalSpace = fopen("./finalPhase.dat", "w");
  FILE *fmag = fopen("./magnet.dat", "w");
  FILE *npout = fopen("./npout.dat", "w");
  FILE *trak = fopen("./trak.dat", "w");
  FILE *virial = fopen("./virial.dat", "w");

  FILE *in = fopen("input.in", "r");
  fscanf(in , "%ld", &n);
  fscanf(in , "%lf", &finalTime);
  fscanf(in , "%lf", &timeStep);
  fscanf(in , "%lf", &timeOutput);
  fscanf(in , "%lf", &p0);
  fscanf(in , "%lf", &r0);
  fscanf(in , "%ld", &seed);
  fclose(in);
  /*
  n = 400000;
  finalTime = 100.;
  timeStep = .1;
  timeOutput = .1;
  p0 = 0.4;
  r0 = pi;
  seed = 10;
  pm = 3.0;
  Np = 32;
  alpha = 2.0;
  idum = -seed;
  */
  /* Inicialização de variáveis*/
  pMean = .0;
  energKin = ran2(&idum);
  energPot = ran2(&idum);
  magX = ran2(&idum);
  magY = ran2(&idum);
  statMoment4 = ran2(&idum);

  dpSize = 2 * pm / Np;
  npMean = .0;
  var = .0;
  Nt = finalTime / timeStep;
  pp = .0;

  double *r = (double *)malloc((double)n * sizeof(double));
  double *p = (double *)malloc((double)n * sizeof(double));
  double *force = (double *)malloc((double)n * sizeof(double));
  long *np = (long *)malloc((long)n * sizeof(long));
  double *state = (double *)malloc((double)n * sizeof(double));
  double *pinit = (double *)malloc((double)n * sizeof(double));
  double *rinit = (double *)malloc((double)n * sizeof(double));
  //  double *gridI = (double *)malloc((double)n * sizeof(double));
  //  double *gridF = (double *)malloc((double)n * sizeof(double));
  // double *tmpP = (double *)malloc((double)n * sizeof(double));
	
  //	double *statMoment4 = (double *)malloc((double)n * sizeof(double));
  //	double *statMoment6 = (double *)malloc((double)n * sizeof(double));

  WaterBag(n, &idum, p0, r0, r, p);
  //	UnitDisk(n, &idum, &massMean, p0, r0, massMin, massMax, r, p, mass);
  //	GaussianPosition(n, &idum, &massMean, p0, r0, massMin, massMax, r, p, mass, alpha);


  for(long i = 0 ; i < n ; ++i)
    {
      pinit[i] = p[i];
      rinit[i] = r[i];
    }
  
  
  /*
    for (long i = 0; i < n; i++)
    {
    np[i] = 0;
    fprintf(init, "%lf\t%lf\n", r[i], p[i]);
    }
  */
  KineticEnergy(n, &energKin, p);
  Force(n, force, r, &magX, &magY);
  PotentialEnergy(n, &energPot, r, magX, magY);
  energ0 = energKin + energPot;

  cout << "Energia Cinetica Inicial: " << energKin << endl;
  cout << "Energia Potencial Inicial: " << energPot << endl;
  cout << "Energia Total Inicial: " << energ0 << endl;
  cout << "Magnetizacoes iniciais:   MagX: " << magX << "  MagY: " << magY << endl;

  error = .0;
  time = .0;
  timeCount = .0;

  long aux = 0;
  double temp = .0f;
  while (time < finalTime)
    {
      temp = .0f;
      /*		
			for (long i = 0; i < n; i++)
			{
			gridI[i] = (long)(p[i] / dpSize);
			}
      */	
      Integration(n, timeStep, &magX, &magY, r, p, force);

      time += timeStep;
      timeCount += timeStep;

      if (timeCount >= timeOutput)
	{
	  KineticEnergy(n, &energKin, p);
	  PotentialEnergy(n, &energPot, r, magX, magY);
	  energ = energKin + energPot;
	  /*		
			for (long i = 0; i < n; i++)
			{
			gridF[i] = (long)(p[i] / dpSize);
			}

			for (long i = 0; i < n; i++) {
			if (gridI[i] != gridF[i]) {
			np[i] += 1;
			}
			npMean += np[i];
			pp += np[i] * np[i];
			}
	  
			statMoment4 = .0;
	  
			for (long i = 0; i < n; i++){
			tmpP[i] += p[i];
			}
	  
			npMean = (npMean/(double)n);
			pp = (pp / ((double)n));

			var = pp - npMean*npMean;
			fprintf(npout, "%lf\t%lf\t%lf\t%lf\n", time, npMean, var, statMoment4);
	  */
	  error = (energ - energ0) / energ0;
	  error = fabs(error);

	  for(long i = 0 ; i < n ; ++i)
	    {
	      temp += p[i]*r[i] - pinit[i]*rinit[i];
	    }
	  temp = 2.*temp/((double) n * time);
	  fprintf(virial,"%lf\t%lf\t%lf\t%lf\n", time, 2*energKin, Virial(n, force, r), temp);

		if((time > 100) && (temp < error)) time = finalTime;
	  
	  //	  printf("%lf\t%1.2le\t%lf\t%lf\t%lf\n", time, error, npMean, var, statMoment4);
	  timeCount = 0.0;
	  //	  fprintf(trak, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", time, p[0], r[0], p[n / 2], r[n / 2], p[n - 1], r[n - 1], p[n/200], r[n/200], p[n / 100], r[n / 100]);
	  fprintf(enrg, "%lf\t%lf\t%lf\n", time, energKin, energPot);
	  fprintf(fmag, "%lf %lf %lf %lf\n", time, magX, magY, sqrt(magX*magX + magY*magY));
	}
    }



  printf("Salvando os espacos de fase finais\n");

  double rr = .0;
  for (long i = 0; i < n; i++)
    {
      rr = r[i];
      while (rr > dpi / 2.){
	rr -= dpi;
      }
      while (rr < -dpi / 2.){
	rr += dpi;
      }

      fprintf(finalSpace, "%lf\t%lf\n", r[i], p[i]);
    }

  free(r);
  free(p);
  free(force);
  free(pinit);
  free(rinit);
  //  free(np);
  //  free(state);
  //  free(gridI);
  //  free(gridF);
  //  free(tmpP);

  fclose(enrg);
  fclose(finalSpace);
  fclose(fmag);
  fclose(init);
  fclose(virial);
  //  fclose(npout);
  //  fclose(trak);
  //	cout << "Gerando graficos" << endl;
  //	system("graph.gp");
  //  system("pause");

  return 0;
}

