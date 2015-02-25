#ifndef HOPPING_H
#define HOPPING_H 1

/* If we  #include <complex.h> before <fftw3.h>, then 
   fftw_complex is the native double complex type. */
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#include "auxiliary.h"
#include "parameter.h"
#include "lattice.h"

double dispersion[L][L][L];

fftw_plan psi1_forward;
fftw_plan psi1_backward;
fftw_plan psi2_forward;
fftw_plan psi2_backward;

void init_fft ();
void hopping (double delta_t);
void clear_fft ();


/***********************************************/
void init_fft ()
{
  psi1_forward = fftw_plan_dft_3d (L, L, L, psi1, 
				   psi1, FFTW_FORWARD, FFTW_MEASURE);
  psi1_backward = fftw_plan_dft_3d (L, L, L, psi1, 
				   psi1, FFTW_BACKWARD, FFTW_MEASURE);
  psi2_forward = fftw_plan_dft_3d (L, L, L, psi2, 
				   psi2, FFTW_FORWARD, FFTW_MEASURE);
  psi2_backward = fftw_plan_dft_3d (L, L, L, psi2, 
				   psi2, FFTW_BACKWARD, FFTW_MEASURE);
  
  for (int i = 0; i < L; ++i)
    for (int j = 0; j < L; ++j)
      for (int k = 0; k < L; ++k)
      {
	if (i <= L/2)
	  dispersion[i][j][k] = -2*t*cos(2*PI*i/L);
	else 
	  dispersion[i][j][k] = -2*t*cos(2*PI*(i-L)/L);
	if (j <= L/2)
	  dispersion[i][j][k] += -2*t*cos(2*PI*j/L);
	else
	  dispersion[i][j][k] += -2*t*cos(2*PI*(j-L)/L);
	if (k <= L/2)
	  dispersion[i][j][k] += -2*t*cos(2*PI*k/L);
	else
	  dispersion[i][j][k] += -2*t*cos(2*PI*(k-L)/L);
      }
}

void hopping (double delta_t)
{
  fftw_execute (psi1_forward);
  fftw_execute (psi2_forward);
  for (int i = 0; i < L; ++i)
    for (int j = 0; j < L; ++j)
      for (int k = 0; k < L; ++k)
      {
	psi1[index(i,j,k)] *= cexp (-I*dispersion[i][j][k]*delta_t)/N;
	psi2[index(i,j,k)] *= cexp (-I*dispersion[i][j][k]*delta_t)/N;
      }
  fftw_execute (psi1_backward);
  fftw_execute (psi2_backward);
}

void clear_fft ()
{
  fftw_destroy_plan (psi1_forward);
  fftw_destroy_plan (psi1_backward);
  fftw_destroy_plan (psi2_forward);
  fftw_destroy_plan (psi2_backward);
}

#endif /* HOPPING_H */
