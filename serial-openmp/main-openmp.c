#include <stdio.h>
#include <stdlib.h>	
#include <math.h>
#include <omp.h>

#include "openmpLibs.h"

#define	pi	3.14159265359
#define dpi	6.28318530718


int main(int argc, char **argv)
{
	long n, seed, idum;
	double p0, r0;
	double energKin, energPot, magX, magY, energ, energ0, error;
	double time, finalTime, timeStep, timeOutput, timeCount;
	

	FILE *init = fopen("./initialPhaseOMP.dat", "w");
	FILE *enrg = fopen("./energyOMP.dat", "w");
	FILE *fmag = fopen("./magnetOMP.dat", "w");
	FILE *finalSpace = fopen("./finalPhaseOMP.dat", "w");

	FILE *in = fopen("./input.in", "r");
	fscanf(in, "%ld", &n);
	fscanf(in, "%lf", &finalTime);
	fscanf(in, "%lf", &timeStep);
	fscanf(in, "%lf", &timeOutput);
	fscanf(in, "%lf", &p0);
	fscanf(in, "%lf", &r0);
	fscanf(in, "%ld", &seed);
	fclose(in);
	
	idum = -seed;

	/* Inicialização de variáveis*/
	energKin = ran2(&idum);
	energPot = ran2(&idum);
	magX = ran2(&idum);
	magY = ran2(&idum);

	double *r = (double *)malloc((double)n * sizeof(double));
	double *p = (double *)malloc((double)n * sizeof(double));
	double *force = (double *)malloc((double)n * sizeof(double));
	
	WaterBag(n, &idum, p0, r0, r, p);

	#pragma omp parallel for
	for (long i = 0; i < n; i++)
	{
		fprintf(init, "%lf\t%lf\n", r[i], p[i]);
	}
	
	KineticEnergy(n, &energKin, p);
	Force(n, force, r, &magX, &magY);
	PotentialEnergy(n, &energPot, r, magX, magY);
	energ0 = energKin + energPot;

	//cout << "Energia Cinetica Inicial: " << energKin << endl;
	//cout << "Energia Potencial Inicial: " << energPot << endl;
	//cout << "Energia Total Inicial: " << energ0 << endl;
	//cout << "Magnetizacoes iniciais:   MagX: " << magX << "  MagY: " << magY << endl;

	error = .0;
	time = .0;
	timeCount = .0;

	while (time < finalTime)
	{
			
		Integration(n, timeStep, &magX, &magY, r, p, force);

		time += timeStep;
		timeCount += timeStep;

		if (timeCount >= timeOutput)
		{
			KineticEnergy(n, &energKin, p);
			PotentialEnergy(n, &energPot, r, magX, magY);
			energ = energKin + energPot;
			error = (energ - energ0) / energ0;
			error = fabs(error);
			// Colocar aqui um if para parar a simulação quandoo errofor grande 
			// Definir erro limite aceitavel
//			printf("%lf\t%1.2le\n", time, error);
			timeCount = 0.0;			
			fprintf(enrg, "%lf\t%lf\t%lf\n", time, energKin, energPot);
			fprintf(fmag, "%lf %lf %lf %lf\n", time, magX, magY, sqrt(magX*magX + magY*magY));
		}
	}

	printf("Salvando os espacos de fase finais\n");

	double rr = .0;
    for (long i = 0; i < n; i++)
	{
		rr = r[i];
		while (rr > dpi / 2.)
		{
			rr -= dpi;
		}
		while (rr < -dpi / 2.)
		{
			rr += dpi;
		}
		fprintf(finalSpace, "%lf\t%lf\n", r[i], p[i]);
	}
 
	free(r);
	free(p);
	free(force);

	fclose(enrg);
	fclose(finalSpace);
	fclose(fmag);
	fclose(init);

	return 0;
}
