#pragma once
#ifdef T8_ENABLE_MRA

namespace t8_mra
{
/// Orthonormal Dubiner scaling functions on the reference triangle (i in
/// [0, 10), tabulated up to P=4). tau1, tau2 are barycentric-ish coords.
double
scaling_function (int i, double tau1, double tau2);

}  // namespace t8_mra
#endif
