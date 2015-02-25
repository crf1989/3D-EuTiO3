#ifndef PARAMETER_H
#define PARAMETER_H 1

/* Units: temperature T in Kelvin, magnetic inductance B in Tesla.
Value of Bohr magneton mu_b is 5.7883818x10^-5 eV/Tesla.  Boltzmman
constant k_B = 8.6173324 x 10^-5 eV /Kelvin.  So we get gnub =
5.788/8.617 = 0.6717, then energy will be measured in Kelvin (more
precisely, Kelvin*k_B is the energy), and magnetic field B in
Tesla. (quote from Prof Wang's program) */

/* parameters for dynamics */
const double dt = 1e-5;
const double t = 0.05/(8.6173324e-5); /* hopping energy */
double JH = 0.05/(8.6173324e-5);

const int samples = 1;
const int time_steps = 20000;
const int stride = 100; 
const int thermalize_steps = 100;


/* the lattice size */
#define L 256
const int N = L*L*L;

/* parameters for lattice */
const double J1 = -0.037;
const double J2 = 0.069;
const double g = 2;
const double u = 0.6717139; /* this is mu_B/k_B */
const double S = 3.968626966596886; /* sqrt (3.5*(3.5+1)) */
const double PI = 3.141592653589793;

double T = 5; /* temperature with default value */
double B = 0; /* magnetic field with default value */


#endif /* PARAMETER_H */
