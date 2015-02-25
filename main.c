#include <stdio.h>
#include <stdlib.h>

#include "parameter.h"
#include "auxiliary.h"
#include "lattice.h"
#include "metropolis.h"
#include "hopping.h"
#include "rotate.h"
#include "trotter_suzuki.h"

int main (int argc, char* argv[])
{
  if (argc != 4)
    {
      perror ("Wrong number of arguments.\n");
      exit (1);
    }
  
  T = atof (argv[1]);
  B = atof (argv[2]);
  JH = atof (argv[3])/(8.6173324e-5);
  printf ("#T = %gK, B = %gT, JH = %geV\n",
	  T, B, atof (argv[3]));

  init ();

  int ts = time_steps/stride + 1;
  double* r2 = calloc (ts, sizeof(double));
  double* r4 = calloc (ts, sizeof(double));
  
  for (int n = 1; n <= samples; ++n)
    {
      thermalize (thermalize_steps); 
      duplicate ();
      
      for (int i = 0; i <= time_steps; ++i)
	{
	  trotter_suzuki_1st (dt);
	  
	  if (i%stride == 0)
	    {
	      double tmp = get_r2 ();
	      r2[i/stride] += tmp;
	      r4[i/stride] += tmp*tmp;
	    }
	}
    }
  
  double deviation;
  for (int i = 0; i < ts; ++i)
    {
      r2[i] /= samples;
      deviation = sqrt ((r4[i] - samples*square(r2[i])) / (samples-1));
      printf ("%g\t%g\t%g\t%g\t%g\n", i*stride*dt, r2[i], deviation, get_norm(),
	      get_spin());
    }
  
  free (r2);
  free (r4);
  clear_all ();

  return 0;      
}
