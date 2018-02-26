#ifndef Virial_h
#define Virial_h
double Virial(long, double *, double *);
double TimeAverage(long, double *, double *);
#endif

double Virial(long n, double *force, double *r)
{
  double aux = .0;
  for(long i = 0 ; i < n ; ++i)
    {
      aux += force[i]*r[i];
    }
  return aux/((double) n);
}


double TimeAvarega(long n, double *p, double *pinit)
{
  double aux = .0;
  for(long i = 0 ; i < n ; ++i)
    {
      aux += p[i]*(p[i] - pinit[i]);
    }

  return aux/((double) n);
}
