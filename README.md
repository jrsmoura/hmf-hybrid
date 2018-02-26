# hmf-hybrid
### Molecular dynamics simulation of a *N*-particles in a cirunference of unitary radius.

Files | Descriotion
------------ | -------------
main-openmp.c | XXXXXX
openmpFuncs.c | XXXXXX
openmpLibs.h  | XXXXXX


**Code Functions**
```C++
void WaterBag(long, long *, double, double, double *, double *);
```
Create a initial phase-space by generating *N* velocities and *N* positions.

```C++
void KineticEnergy(long, double *, double *);
```
Calculates the kinetic energy per particle from the system.

```C++
void PotentialEnergy(long, double *, double *, double, double);
```
Calculates the potential energy per particle from the system.

```C++
void Force(long, double *, double *, double *, double *);
```


```
void Integration(long, double, double *, double *, double *, double *, double *);
```

```C++
float ran2(long *);
```
Generate a uniform distributed random number between 0 and 1.
