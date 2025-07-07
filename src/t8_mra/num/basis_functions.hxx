#pragma once
#ifdef T8_ENABLE_MRA

namespace t8_mra
{
double
skalierungsfunktion (int i, double tau1, double tau2);

double
muttermultiwavelets (int p, int i, double tau1, double tau2, int e);

double
muttermultiwavelet (int p, int i, double tau1, double tau2);

double
skalierungsfunktion_nextlevel (int i, double tau1, double tau2);

}  // namespace t8_mra
#endif
