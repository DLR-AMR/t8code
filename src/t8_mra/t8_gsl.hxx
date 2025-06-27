#ifndef T8_GSL_HXX
#define T8_GSL_HXX

#include <t8.h> /* General t8code header, always include this. */
#include <iostream>
#include <fstream>
#ifdef T8_ENABLE_MRA
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

typedef double (*func) (double, double);
typedef double (*spline) (const gsl_spline2d* spline, gsl_interp_accel* xacc, gsl_interp_accel* yacc);

// Function to read MPTRAC data
void
read_mptrac_data (const std::string& filename, float* fileData, size_t dataSize);

// Function to initialize the spline
void
setup_spline (const float* fileData, size_t nx, size_t ny, gsl_spline2d* spline, gsl_interp_accel* xacc,
              gsl_interp_accel* yacc, double* za);

// Function to allocate memory for spline resources
void
allocate_spline_resources (const gsl_interp2d_type* T, size_t nx, size_t ny, gsl_spline2d** spline,
                           gsl_interp_accel** xacc, gsl_interp_accel** yacc, double** za);

// Function to clean up allocated memory
void
free_spline_resources (gsl_spline2d* spline, gsl_interp_accel* xacc, gsl_interp_accel* yacc, double* za);

// Function to evaluate the spline
double
EvaluateSpline (const gsl_spline2d* spline, const double x, const double y, gsl_interp_accel* xacc,
                gsl_interp_accel* yacc);

#endif
#endif
