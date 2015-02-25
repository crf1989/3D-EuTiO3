#ifndef ROTATE_H
#define ROTATE_H 1

#include <complex.h>
#include <math.h>

#include "lattice.h"

/* rotate_electron include the effect of Zeeman term and spin term */
void rotate_electron (double delta_t);

/* rotate_spin include the effect of Zeeman term, electron term and
   spin term */
void rotate_spin_sublattice (int xbegin, int ybgein, int zbegin,
			     int xinc, int yinc, int zinc, double delta_t);


/*******************************************************/
void rotate_electron (double delta_t)
/* This function rotate the spinor, the formula of exponential of
   Pauli vector is used, see
   http://en.wikipedia.org/wiki/Pauli_matrices */
{ 
  double complex sigma[2][2];
  
  for (int i = 0; i < N; ++i)
    {
      double s[3] = {0, 0, 0};
      for (int j = 0; j < 8; ++j)
	for (int k = 0; k < 3; ++k)
	  s[k] += 0.5*JH*copy[elatt[i].sNN[j]].s[k];
      s[2] += u*B;

      double phi = sqrt (square(s[0]) + square(s[1]) + square(s[2]));
      
      /* if |s| is zero, do nothing. */
      if (fabs (phi) < 1e-8)
	return;

      /* \vec{s}.\vec{simga}/|s|, here simga is Pauli matrice. */
      for (int k = 0; k < 3; ++k)
	s[k] /= phi;

      double cth = cos (phi*delta_t);
      double sth = sin (phi*delta_t);
      
      sigma[0][0] = cth + I*s[2]*sth;
      sigma[0][1] = I*s[0]*sth + s[1]*sth;
      sigma[1][0] = I*s[0]*sth - s[1]*sth;
      sigma[1][1] = cth - I*s[2]*sth;

      /* matrice multiplication */
      double complex tmp1 = psi1[i];
      double complex tmp2 = psi2[i];
      
      psi1[i] = sigma[0][0]*tmp1 + sigma[0][1]*tmp2;
      psi2[i] = sigma[1][0]*tmp1 + sigma[1][1]*tmp2;
    }
}


void rotate_spin_sublattice (int xbegin, int ybegin, int zbegin,
			     int xinc, int yinc, int zinc, double delta_t)
/* The 3 dimensional rotation matrix see 
   http://en.wikipedia.org/wiki/Rotation_matrix */
{
  double R[3][3];
  for (int i = xbegin; i < L && i >= 0; i += xinc)
    for (int j = ybegin; j < L && j >= 0; j += yinc)
      for (int k = zbegin; k < L && k >= 0; k += zinc)
      {
	double H[3] = {0, 0, 0};
	int m = index (i,j,k);
	/* nearest neighbors */
	for (int n = 0; n < 6; ++n)
	  {
	    H[0] += J1 * copy[slatt[m].sNN[n]].s[0];
	    H[1] += J1 * copy[slatt[m].sNN[n]].s[1];
	    H[2] += J1 * copy[slatt[m].sNN[n]].s[2];
	  }
	/* next nearest neighbors */
	for (int n = 0; n < 12; ++n)
	  {
	    H[0] += J2 * copy[slatt[m].sNNN[n]].s[0];
	    H[1] += J2 * copy[slatt[m].sNNN[n]].s[1];
	    H[2] += J2 * copy[slatt[m].sNNN[n]].s[2];
	  }
	/* Zeeman term */
	H[2] += g*u*B;
	/* eletron term */
	for (int n = 0; n < 8; ++n)
	  {
	    H[0] += JH * creal (conj(psi1[slatt[m].eNN[n]])*psi2[slatt[m].eNN[n]]);
	    H[1] += JH * cimag (conj(psi1[slatt[m].eNN[n]])*psi2[slatt[m].eNN[n]]);
	    H[2] += JH * 0.5*(cabs2 (psi1[slatt[m].eNN[n]]) - 
			  cabs2 (psi2[slatt[m].eNN[n]]));
	  }
	/* angular velocity */
	double omega = sqrt (square(H[0]) + square(H[1]) + square(H[2]));
	/* unit rotation axis */
	if (fabs (omega) < 1e-8)
	  continue;
	H[0] /= omega; H[1] /= omega; H[2] /= omega;
	double cth = cos (-delta_t*omega); /* cos (theta) */
	double sth = sin (-delta_t*omega); /* sin (theta) */

	R[0][0] = cth + H[0]*H[0]*(1-cth);
	R[0][1] = H[0]*H[1]*(1-cth) - H[2]*sth;
	R[0][2] = H[0]*H[2]*(1-cth) + H[1]*sth;
	
	R[1][0] = H[1]*H[0]*(1-cth) + H[2]*sth;
	R[1][1] = cth + H[1]*H[1]*(1-cth);
	R[1][2] = H[1]*H[2]*(1-cth) - H[0]*sth;
	
	R[2][0] = H[2]*H[0]*(1-cth) - H[1]*sth;
	R[2][1] = H[2]*H[1]*(1-cth) + H[0]*sth;
	R[2][2] = cth + H[2]*H[2]*(1-cth);

	double stx = copy[m].s[0];
	double sty = copy[m].s[1];
	double stz = copy[m].s[2];
	
	/* matrix multiplication */
	for (int k = 0; k < 3; ++k)
	  copy[m].s[k] = R[k][0]*stx + R[k][1]*sty + R[k][2]*stz;
      }
}

void rotate_spin_forward (double delta_t)
{
  rotate_spin_sublattice (0, 0, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 1, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 0, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 1, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 0, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 1, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 0, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 1, 1, 2, 2, 2, delta_t);
}

void rotate_spin_backward (double delta_t)
{
  rotate_spin_sublattice (1, 1, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 0, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 1, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 0, 1, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 1, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (1, 0, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 1, 0, 2, 2, 2, delta_t);
  rotate_spin_sublattice (0, 0, 0, 2, 2, 2, delta_t);
}

/* void rotate_spin_backward (double delta_t) */
/* { */
/*   rotate_spin_sublattice (L-1, L-1, L-1, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-1, L-2, L-1, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-2, L-1, L-1, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-2, L-2, L-1, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-1, L-1, L-2, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-1, L-2, L-2, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-2, L-1, L-2, -2, -2, -2, delta_t); */
/*   rotate_spin_sublattice (L-2, L-2, L-2, -2, -2, -2, delta_t); */
/* } */

#endif /* ROTATE_H */
