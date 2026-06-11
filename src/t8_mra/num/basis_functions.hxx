#pragma once
#ifdef T8_ENABLE_MRA

namespace t8_mra
{
double
scaling_function (int i, double tau1, double tau2);

double
muttermultiwavelets (int p, int i, double tau1, double tau2, int e);

double
muttermultiwavelet (int p, int i, double tau1, double tau2);

double
scaling_function_nextlevel (int i, double tau1, double tau2);

}  // namespace t8_mra
#endif
