#ifndef __BASIS_FUNCTIONS_INCLUDE__
#define __BASIS_FUNCTIONS_INCLUDE__

#include <cmath>

double
skalierungsfunktion (int i, double tau1, double tau2);

double
muttermultiwavelets (int p, int i, double tau1, double tau2, int e);

double
muttermultiwavelet (int p, int i, double tau1, double tau2);

double
skalierungsfunktion_nextlevel (int i, double tau1, double tau2);
#endif
