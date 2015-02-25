#ifndef TROTTER_SUZUKI_H
#define TROTTER_SUZUKI_H 1

#include "hopping.h"
#include "rotate.h"

/* this init function initialize all things */
void init ();
void clear_all ();

void trotter_suzuki_1st (double delta_t);
void trotter_suzuki_2nd (double delta_t);
void trotter_suzuki_4th (double delta_t);


/******************************************************/
void init ()
{
  init_lattice ();
  init_workspace ();
  init_fft ();
  init_metropolis ();
}

void clear_all ()
{
  clear_lattice ();
  clear_workspace ();
  clear_fft ();
  clear_metroplis ();
}

void trotter_suzuki_1st (double delta_t)
{
  rotate_spin_forward (delta_t);
  rotate_electron (delta_t);
  hopping (delta_t);
}

void trotter_suzuki_2nd (double delta_t)
{
  rotate_spin_forward (delta_t/2);
  rotate_electron (delta_t/2);
  hopping (delta_t);
  rotate_electron (delta_t/2);
  rotate_spin_backward (delta_t/2);
}

void trotter_suzuki_4th (double delta_t)
{
  trotter_suzuki_2nd (1.351207191959657*delta_t);
  trotter_suzuki_2nd ((1-2*1.351207191959657)*delta_t);
  trotter_suzuki_2nd (1.351207191959657*delta_t);
}

#endif /* TROTTER_SUZUKI_H */
