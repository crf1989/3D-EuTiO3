#ifndef METROPOLIS_H
#define METROPOLIS_H 1

#include <stdlib.h>
#include <time.h>
#include "auxiliary.h"
#include "parameter.h"
#include "lattice.h"

void init_metropolis ();
void metropolis (int steps);
void thermalize (int steps);
void clear_metroplis ();


/**********************************************/
void init_metropolis ()
{
  srand (time (0));
}

void metropolis (int steps)
{
  int m;
  double sn[3]; /* new spin vector */
  double de; /* delta energy */

  while (steps--)
    {
      m = rand () % N;
      genvec (&sn[0], &sn[1], &sn[2]);
      
      /* energy change due to magnetic field */
      de = -g*u*B * (sn[2] - slatt[m].s[2]);
      
      /* energy change due to nearest neighbors */
      for (int i = 0; i < 6; ++i)
	for (int j = 0; j < 3; ++j)
	  de += -J1 * (sn[j] - slatt[m].s[j]) * slatt[slatt[m].sNN[i]].s[j];
      
      /* energy change due to next nearest neighbors */
      for (int i = 0; i < 12; ++i)
	for (int j = 0; j < 3; ++j)
	  de += -J2 * (sn[j] - slatt[m].s[j]) * slatt[slatt[m].sNNN[i]].s[j];

      if (de <= 0 || drand() < exp (-de/T))
	for (int i = 0; i < 3; ++i)
	  slatt[m].s[i] = sn[i];
    }
}

void thermalize (int steps)
{
  while (steps--)
    metropolis (N);
}

void clear_metroplis ()
{
  return;
}


#endif /* METROPOLIS_H */
