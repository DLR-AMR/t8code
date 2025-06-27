#include <iostream>
#include <fstream>
#include <t8.h>

#ifdef T8_ENABLE_MRA
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

void
read_mptrac_data (const std::string& filename, float* fileData, size_t dataSize)
{
  std::ifstream inputFileStream (filename, std::ios::in | std::ios::binary);
  if (!inputFileStream) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  std::cout << "Reading binary data from: " << filename << std::endl;
  inputFileStream.read ((char*) fileData, dataSize * sizeof (float));
  inputFileStream.close ();
}

// Function to initialize the spline
void
setup_spline (const float* fileData, size_t nx, size_t ny, gsl_spline2d* spline, gsl_interp_accel* xacc,
              gsl_interp_accel* yacc, double* za)
{
  double xa[nx], ya[ny];

  // Define grid points
  for (size_t i = 0; i < nx; ++i) {
    xa[i] = 0.3 * i;
  }
  for (size_t i = 0; i < ny; ++i) {
    ya[i] = -90.0 + 0.3 * i;
  }

  // Set the z grid values from fileData
  for (size_t iy = 600, i = 0; iy >= 0 && i < ny; --iy, ++i) {
    for (size_t ix = 0; ix < nx; ++ix) {
      gsl_spline2d_set (spline, za, ix, i, fileData[ix * ny + iy]);
    }
  }

  // Initialize the spline object
  gsl_spline2d_init (spline, xa, ya, za, nx, ny);
}

// Function to allocate memory for spline resources
void
allocate_spline_resources (const gsl_interp2d_type* T, size_t nx, size_t ny, gsl_spline2d** spline,
                           gsl_interp_accel** xacc, gsl_interp_accel** yacc, double** za)
{
  // Allocate memory for the z array (za) using T8_ALLOC
  *za = T8_ALLOC (double, nx* ny);

  // Allocate memory for the spline object using gsl_spline2d_alloc
  *spline = gsl_spline2d_alloc (T, nx, ny);

  // Allocate memory for the x and y accelerators using gsl_interp_accel_alloc
  *xacc = gsl_interp_accel_alloc ();
  *yacc = gsl_interp_accel_alloc ();
}

// Function to clean up allocated memory
void
free_spline_resources (gsl_spline2d* spline, gsl_interp_accel* xacc, gsl_interp_accel* yacc, double* za)
{
  gsl_spline2d_free (spline);    // Free the spline object
  gsl_interp_accel_free (xacc);  // Free x accelerator
  gsl_interp_accel_free (yacc);  // Free y accelerator
  T8_FREE (za);                  // Free the dynamically allocated z array (za)
}

double
EvaluateSpline (const gsl_spline2d* spline, const double x, const double y, gsl_interp_accel* xacc,
                gsl_interp_accel* yacc)
{
  return gsl_spline2d_eval (spline, x, y, xacc, yacc);
}
#endif
