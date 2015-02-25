#ifndef LATTICE_H
#define LATTICE_H 1

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "parameter.h"
#include "auxiliary.h"

typedef struct 
{
  double s[3];
  int sNN[6];
  int sNNN[12];
  int eNN[8];
} ssite;

typedef struct 
{
  int eNN[6];
  int sNN[8];
} esite;

ssite* slatt; /* spin lattice */
esite* elatt; /* electron lattice */

int index (int i, int j, int k);
void position (int m, int* i, int* j, int* k);
void init_lattice ();
void clear_lattice ();

/* workspace for dynamics */
double complex* psi1; /* spin up */
double complex* psi2; /* spin down */

typedef struct
{
  double s[3];
} workspace;
workspace* copy;

void init_workspace ();
void duplicate ();
void clear_workspace ();

/* functions to check */
double get_r2 ();
double get_norm ();
double get_spin ();

int index (int i, int j, int k)
{
  if (i < 0) i += L;
  if (i >= L) i -= L;
  if (j < 0) j += L;
  if (j >= L) j -= L;
  if (k < 0) k += L;
  if (k >= L) k -= L;
  return i*L*L + j*L + k;
}

void position (int m, int* i, int* j, int* k)
{
  *k = m % L;
  *j = (m % (L*L)) / L;
  *i = m / (L*L);
}

/* initialize the lattice and neighbors */
void init_lattice ()
{
  slatt = malloc (N * sizeof(ssite)); 
  elatt = malloc (N * sizeof(esite));
  
  for (int i = 0; i < N; ++i)
    {
      /* initialize the spin */
      genvec (&slatt[i].s[0], &slatt[i].s[1], &slatt[i].s[2]);

      int a, b, c;
      position (i, &a, &b, &c);
      /* set the nearest neighbors between spin and spin */
      slatt[i].sNN[0] = index (a-1, b, c);
      slatt[i].sNN[1] = index (a+1, b, c);
      slatt[i].sNN[2] = index (a, b-1, c);
      slatt[i].sNN[3] = index (a, b+1, c);
      slatt[i].sNN[4] = index (a, b, c-1);
      slatt[i].sNN[5] = index (a, b, c+1);
      /* set the next nearest neighbors between spin and spin */
      slatt[i].sNNN[0] = index (a-1, b-1, c);
      slatt[i].sNNN[1] = index (a-1, b+1, c);
      slatt[i].sNNN[2] = index (a+1, b-1, c);
      slatt[i].sNNN[3] = index (a+1, b+1, c);
      slatt[i].sNNN[4] = index (a-1, b, c-1);
      slatt[i].sNNN[5] = index (a-1, b, c+1);
      slatt[i].sNNN[6] = index (a+1, b, c-1);
      slatt[i].sNNN[7] = index (a+1, b, c+1);
      slatt[i].sNNN[8] = index (a, b-1, c-1);
      slatt[i].sNNN[9] = index (a, b-1, c+1);
      slatt[i].sNNN[10] = index (a, b+1, c-1);
      slatt[i].sNNN[11] = index (a, b+1, c+1);
      /* set the nearest neighbors between electron and electron */
      elatt[i].eNN[0] = index (a-1, b, c);
      elatt[i].eNN[1] = index (a+1, b, c);
      elatt[i].eNN[2] = index (a, b-1, c);
      elatt[i].eNN[3] = index (a, b+1, c);
      elatt[i].eNN[4] = index (a, b, c-1);
      elatt[i].eNN[5] = index (a, b, c+1);
      
      /* set the nearest neighbors between electron and spin */
      elatt[i].sNN[0] = index (a, b, c);
      elatt[i].sNN[1] = index (a+1, b, c);
      elatt[i].sNN[2] = index (a, b-1, c);
      elatt[i].sNN[3] = index (a+1, b-1, c);
      elatt[i].sNN[4] = index (a, b, c-1);
      elatt[i].sNN[5] = index (a+1, b, c-1);
      elatt[i].sNN[6] = index (a, b-1, c-1);
      elatt[i].sNN[7] = index (a+1, b-1, c-1);

      slatt[i].eNN[0] = index (a, b, c);
      slatt[i].eNN[1] = index (a-1, b, c);
      slatt[i].eNN[2] = index (a, b+1, c);
      slatt[i].eNN[3] = index (a-1, b+1, c);
      slatt[i].eNN[4] = index (a, b, c+1);
      slatt[i].eNN[5] = index (a-1, b, c+1);
      slatt[i].eNN[6] = index (a, b+1, c+1);
      slatt[i].eNN[7] = index (a-1, b+1, c+1);
    }
}

void clear_lattice ()
{
  free (slatt);
  free (elatt);
}

void init_workspace ()
{
  psi1 = malloc (N * sizeof(double complex));
  psi2 = malloc (N * sizeof(double complex));
  copy = malloc (N * sizeof(workspace));
}

void duplicate ()
{
  for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < 3; ++j)
	copy[i].s[j] = slatt[i].s[j];
      psi1[i] = 0;
      psi2[i] = 0;
    }
  psi1[index(L/2,L/2,L/2)] = 1;
}

void clear_workspace ()
{
  free (psi1);
  free (psi2);
  free (copy);
}

double get_r2 ()
{
  double result = 0;
  for (int i = 0; i < L; ++i)
    for (int j = 0; j < L; ++j)
      for (int k = 0; k < L; ++k)
	result += (square(i-L/2) + square(j-L/2) + square(k-L/2)) * 
	  (cabs2 (psi1[index(i,j,k)]) + cabs2 (psi2[index(i,j,k)]));

  return result;
}

double get_norm ()
{
  double result = 0;
  for (int i = 0; i < N; ++i)
    result += cabs2 (psi1[i]) + cabs2 (psi2[i]);

  return result;
}

double get_spin ()
{
  double result = 0;
  for (int i = 0; i < N; ++i)
    result += sqrt (square (copy[i].s[0]) +
		    square (copy[i].s[1]) + 
		    square (copy[i].s[2]));
  
  return result/N;
}

#endif /* LITTICE_H */
