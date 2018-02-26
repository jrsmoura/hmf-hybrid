# Hamiltonian Mean Field simulation using OpenMP

### Molecular dynamics simulation of a *N*-particles in a cirunference of unitary radius.

Files | Descriotion
------------ | -------------
main-openmp.c | Main part of the code
openmpFuncs.c | Contains all functions
openmpLibs.h  | Functions headers


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
Calculates the mean force from the system

```C++
void Integration(long, double, double *, double *, double *, double *, double *);
```
Integrate equations of motion using a sympletic integrator, fourth order Yoshida.

```C++
float ran2(long *);
```
Generate a uniform distributed random number between 0 and 1 (Numerical Recipies C/C++, 3rd Ed.)


## The figures below show some basics results that can be obtain

![Time evolution of inetic and potential energies per particle](/hmf-hybrid/Energies.png)
Format: ![Alt Text](url)
