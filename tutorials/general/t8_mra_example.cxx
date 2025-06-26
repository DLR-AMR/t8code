#include <t8.h>                          /* General t8code header, always include this. */
#include <t8_cmesh.hxx>                  /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h> /* forest definition and basic interface. */
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_io.h>                                       /* forest io interface. */
#include <t8_schemes/t8_default/t8_default.hxx>                           /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>                              /* geometrical information */
#include "t8_mra/vecmat.hxx"
#include "t8_mra/basis_functions.hxx"
#include "t8_mra/mask_coefficients.hxx"
#include "t8_mra/dunavant.hxx"
#include "t8_mra/t8_gsl.hxx"
#include "t8_mra/t8_mra.hxx"
#include "t8_mra/t8_unordered_dense.hxx"
#include <iomanip>
// #include "t8_unordered_dense.hxx"
#include <cmath>
#include <vector>
#include <sc_statistics.h>
#include <t8_refcount.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_trees.h>
//#include <t8_element_c_interface.h>
#include <iostream>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_interp2d.h>
// #include <gsl/gsl_spline2d.h>
#include <fstream>
#include <time.h>
// 3. Declare function prototypes or class declarations (if needed)

// 4. Optionally, include any namespaces you want to use (avoid if using standard practice)
using namespace std;

/* We build 4 kinds of cmeshes in the following by hand. */

/* This cmesh is [0,1]² with 8 triangles. */
t8_cmesh_t
t8_cmesh_new_basic (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[72] = {
    0,   0,   0, 0.5, 0,   0, 0.5, 0.5, 0,  //triangle 2
    0,   0,   0, 0,   0.5, 0, 0.5, 0.5, 0,  //triangle 1
    0.5, 0,   0, 1,   0,   0, 1,   0.5, 0,  //triangle 3
    0.5, 0,   0, 0.5, 0.5, 0, 1,   0.5, 0,  //triangle 4
    0,   0.5, 0, 0.5, 0.5, 0, 0.5, 1,   0,  //triangle 5
    0,   0.5, 0, 0,   1,   0, 0.5, 1,   0,  //triangle 6
    0.5, 0.5, 0, 1,   0.5, 0, 1,   1,   0,  //triangle 7
    0.5, 0.5, 0, 0.5, 1,   0, 1,   1,   0,  //triangle 8
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 4, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 6, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 7, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 1, 1, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is [0,1]² with 2 triangles. for debugging purposes to compare to results from Florian Sieglar */
t8_cmesh_t
t8_cmesh_new_debugging (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[18] = {
    0, 0, 0, 0, 1, 0, 1, 0, 0,  //triangle 1
    1, 1, 0, 1, 0, 0, 0, 1, 0,  //triangle 2
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is an octagon with 14 elements */
t8_cmesh_t
t8_cmesh_new_octagon (sc_MPI_Comm comm)
{
  double a = 2.4142135623731;  //1+sqrt(2)
  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[126] = {
    -1, -a, 0, -1, -1, 0, -a, -1, 0,  //triangle 1
    -1, -a, 0, 1,  -a, 0, 1,  -1, 0,  //triangle 2
    -1, -a, 0, -1, -1, 0, 1,  -1, 0,  //triangle 3
    1,  -a, 0, 1,  -1, 0, a,  -1, 0,  //triangle 4
    -a, -1, 0, -1, -1, 0, -1, 1,  0,  //triangle 5
    -a, -1, 0, -a, 1,  0, -1, 1,  0,  //triangle 6
    -1, -1, 0, 1,  -1, 0, 1,  1,  0,  //triangle 7
    -1, -1, 0, -1, 1,  0, 1,  1,  0,  //triangle 8
    1,  -1, 0, a,  -1, 0, 1,  1,  0,  //triangle 9
    a,  -1, 0, 1,  1,  0, a,  1,  0,  //triangle 10
    -a, 1,  0, -1, 1,  0, -1, a,  0,  //triangle 11
    -1, 1,  0, 1,  a,  0, 1,  1,  0,  //triangle 12
    -1, 1,  0, -1, a,  0, 1,  a,  0,  //triangle 13
    1,  a,  0, 1,  1,  0, a,  1,  0,  //triangle 14
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 8, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 9, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 10, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 11, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 12, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 13, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);
  t8_cmesh_set_tree_vertices (cmesh, 8, vertices + 72, 3);
  t8_cmesh_set_tree_vertices (cmesh, 9, vertices + 81, 3);
  t8_cmesh_set_tree_vertices (cmesh, 10, vertices + 90, 3);
  t8_cmesh_set_tree_vertices (cmesh, 11, vertices + 99, 3);
  t8_cmesh_set_tree_vertices (cmesh, 12, vertices + 108, 3);
  t8_cmesh_set_tree_vertices (cmesh, 13, vertices + 117, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 2, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 4, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 6, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 8, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 7, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 10, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 8, 0, 1, 1);
  t8_cmesh_set_join (cmesh, 7, 11, 0, 1, 1);
  t8_cmesh_set_join (cmesh, 8, 9, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 9, 13, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 10, 12, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 11, 12, 2, 1, 1);
  t8_cmesh_set_join (cmesh, 11, 13, 0, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is a complex polygonal shape with 8 triangles. */
t8_cmesh_t
t8_cmesh_new_complex_polygonal_shape (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[72] = {
    0.4, 0.5, 0, 1,   0,   0, 1,   0.6, 0,  //triangle 1
    0.4, 0.2, 0, 0.4, 0.5, 0, 1,   0,   0,  //triangle 2
    0,   0,   0, 0.4, 0.2, 0, 0.4, 0.5, 0,  //triangle 3
    0,   0,   0, 0.2, 0.4, 0, 0.4, 0.5, 0,  //triangle 4
    0.4, 0.8, 0, 0.2, 0.4, 0, 0.4, 0.5, 0,  //triangle 5
    0.4, 0.8, 0, 0.2, 0.4, 0, 0,   0.6, 0,  //triangle 6
    0.4, 0.8, 0, 0.4, 1,   0, 0,   0.6, 0,  //triangle 7
    0.4, 0.8, 0, 0.4, 1,   0, 1,   0.8, 0,  //triangle 8
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 2, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 6, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 2, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is an L-shape with 4 triangles. */
t8_cmesh_t
t8_cmesh_new_l_shape (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[36] = {
    0.5, 0.5, 0, 1, 0, 0, 1,   0.5, 0,  //triangle 1
    0.5, 0.5, 0, 1, 0, 0, 0,   0,   0,  //triangle 2
    0.5, 0.5, 0, 0, 1, 0, 0,   0,   0,  //triangle 3
    0.5, 0.5, 0, 0, 1, 0, 0.5, 1,   0,  //triangle 4
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 2, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is (x-axis)[0,360]x(y-axis)[-90,90] with 4 triangles (specifically for the era5 data). */
t8_cmesh_t
t8_cmesh_new_earth (sc_MPI_Comm comm)
{

  double vertices[36] = {
    0,   -90, 0, 180, -90, 0, 180, 90,  0,  //triangle 1
    360, 90,  0, 180, -90, 0, 360, -90, 0,  //triangle 2
    0,   -90, 0, 0,   90,  0, 180, 90,  0,  //triangle 3
    360, 90,  0, 180, -90, 0, 180, 90,  0   //triangle 4
  };

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);

  t8_cmesh_set_join (cmesh, 0, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 2, 2, 0);

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* These are six test functions. */
double
f1 (double x, double y)
{
  //if ((x == -1.) && (y == -1.)) return 1.;
  return sin (2. * M_PI * x) * sin (2. * M_PI * y);
}

double
f2 (double x, double y)
{
  //if ((x == -1.) && (y == -1.)) return 2.;
  double r = x * x + y * y;
  return (r < 0.25) ? 1.0 : 0.0;
}

double
f3 (double x, double y)
{
  //if ((x == -1.) && (y == -1.)) return 3.;
  double r = x * x + y * y;
  return (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x);
}

double
f5 (double x, double y)
{
  //if ((x == -1.) && (y == -1.)) return 4.;
  return sin (1 / (1.001 - x * y));
}

double
f4 (double x, double y)
{
  //if ((x == -1.) && (y == -1.)) return 6.;
  if (x < 0.41)
    return 0.;
  double r4 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
  double r = sqrt (r4);
  if (r > 1. / 3.)
    return 0.;
  r *= 3.;
  r4 *= 9.;
  r4 *= r4;
  double rm1 = r - 1.;
  double rm1h2 = rm1 * rm1;
  double rm1h3 = rm1 * rm1h2;
  return 1. - r4 + 4. * r4 * rm1 - 10. * r4 * rm1h2 + 20 * r4 * rm1h3;
}

/* Testing jump sizes */
double
f6 (double x, double y)
{
  double val = 0.0;
  if ((y >= 0.1))
    val += 1e-8;
  if ((y >= 0.2))
    val += 1e-7;
  if ((y >= 0.3))
    val += 1e-6;
  if ((y >= 0.4))
    val += 1e-5;
  if ((y >= 0.5))
    val += 1e-4;
  if ((y >= 0.6))
    val += 1e-3;
  if ((y >= 0.7))
    val += 1e-2;
  if ((y >= 0.8))
    val += 1e-1;
  if ((y >= 0.9))
    val += 1.0;
  return val;
}

double
f7 (double x, double y)
{
  return y * y * y * y * y * y * y * x * x * x * x * x * x * x + x * y - x + 10;
}

//The polynomials are used to verify the Vanishing Moments
double
zero_degree_const (double x, double y)
{
  return 1;
}

//The polynomials are used to verify the Vanishing Moments
double
zero_degree_large_const (double x, double y)
{
  return 100000000000;
}

double
first_degree_x (double x, double y)
{
  return x;
}

double
first_degree_y (double x, double y)
{
  return y;
}

double
second_degree_x (double x, double y)
{
  return x * x;
}

double
second_degree_y (double x, double y)
{
  return y * y;
}

double
second_degree_xy (double x, double y)
{
  return x * y;
}

double
third_degree_x (double x, double y)
{
  return x * x * x;
}

double
third_degree_xxy (double x, double y)
{
  return x * x * y;
}

double
third_degree_xyy (double x, double y)
{
  return x * y * y;
}

double
third_degree_y (double x, double y)
{
  return y * y * y;
}

double
rho (double x, double y)
{
  // Constants
  double x0 = 0.5, y0 = 0.5;
  double dx = x - x0;
  double dy = y - y0;

  double r = sqrt (dx * dx + dy * dy);
  double theta = atan2 (dy, dx);

  // Core Gaussian jet
  double core = exp (-80.0 * (dx * dx + dy * dy));

  // Radial ripples
  double ripples = sin (30.0 * r - 5.0 * theta) * exp (-40.0 * r);

  // Shock waves
  double shock_x = sin (20.0 * M_PI * x) * exp (-100.0 * (y - 0.5) * (y - 0.5));
  double shock_y = sin (20.0 * M_PI * y) * exp (-100.0 * (x - 0.5) * (x - 0.5));

  // Diagonal interaction terms
  double diagonal1 = exp (-300.0 * (x - y - 0.1) * (x - y - 0.1));
  double diagonal2 = exp (-300.0 * (x + y - 1.1) * (x + y - 1.1));

  // Swirl term
  double swirl = sin (10.0 * theta) * exp (-20.0 * r);

  // Combine all components
  double result
    = 0.5 * core + 0.3 * ripples + 0.2 * shock_x + 0.2 * shock_y + 0.15 * diagonal1 + 0.15 * diagonal2 + 0.2 * swirl;

  return result;
}

double
complex_flow (double x, double y)
{
  // Center and shift coordinates
  double cx = 0.5, cy = 0.5;
  double dx = x - cx;
  double dy = y - cy;
  double r = sqrt (dx * dx + dy * dy);
  double theta = atan2 (dy, dx);

  // Turbulence-like small scale noise (pseudo-turbulence)
  double turbulence = sin (70 * x) * cos (70 * y) * 0.1;

  // Ripple pattern (damped radial sine)
  double ripples = sin (25.0 * r - 4.0 * theta) * exp (-30.0 * r);

  // Spiral vortex (decaying angular swirl)
  double vortex = sin (12.0 * theta + 8.0 * r) * exp (-10.0 * r);

  // Jump/shock-like layer (sharp sigmoidal transitions)
  double jump_x = 1.0 / (1.0 + exp (-100.0 * (x - 0.5)));
  double jump_y = 1.0 / (1.0 + exp (-100.0 * (y - 0.7)));
  double jumps = jump_x * (1.0 - jump_y);

  // Shear layer with sine disturbance
  double shear = tanh (20.0 * (x - y)) * sin (10.0 * (x + y));

  // Combine all effects
  double result = 0.3 * turbulence + 0.3 * ripples + 0.2 * vortex + 0.3 * jumps + 0.2 * shear;

  return result;
}

double
fractal_field (double x, double y)
{
  // Map [0,1]^2 to [-2,2]^2 for typical fractal behavior
  double zx = 4.0 * x - 2.0;
  double zy = 4.0 * y - 2.0;

  double cx = zx;
  double cy = zy;

  int max_iter = 1000;
  int iter = 0;

  while (zx * zx + zy * zy < 4.0 && iter < max_iter) {
    double xtemp = zx * zx - zy * zy + cx;
    zy = 2.0 * zx * zy + cy;
    zx = xtemp;
    iter++;
  }

  // Normalize to [0,1] range
  return (double) iter / max_iter;
}

double
testtest (double x, double y)
{
  // Avoid division by zero or tan(y) singularities
  double tan_y = tan (y);
  if (fabs (tan_y) < 1e-8 || fabs (1.0 - y * y) < 1e-8) {
    return 0.0;  // or use NAN / DBL_MAX to signal undefined
  }

  double exp_part = exp (3.0 * x * x);
  double sin_part = sin (exp_part / tan_y);
  double cos_part = cos (1.0 / (1.0 - y * y));
  double final = sin_part * cos_part * exp (y);

  return final;
}

// /** Write the forest and the data corresponding to the forest in a vtu file.
//  *
//  * \param [in] forest           Forest that should written in the vtu file
//  * \param [in] data             Data corresponding to the forest
//  * \param [in] prefix           Prefix to define the file name
//  */
// static void
// t8_write_vtu_mra (t8_forest_t forest, const char *prefix)
// {
//   const t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
//   /* We need to allocate a new array to store the volumes on their own.
//    * This array has one entry per local element. */
//   double *element_data = T8_ALLOC (double, num_elements);
//   /* The number of user-defined data fields to write. */
//   int num_data = 1;
//   t8_vtk_data_field_t vtk_data;
//   vtk_data.type = T8_VTK_SCALAR;
//   strcpy (vtk_data.description, "Element own data");
//   vtk_data.data = element_data;
//   /* Copy the element's data from the data array to the output array. */
//   for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
//     element_data[ielem] = t8_element_get_value (data, ielem).values;
//   }
//
//   /* To write user-defined data, we need to extend the output function t8_forest_vtk_write_file
//    * from t8_forest_vtk.h. Despite writing user data, it also offers more control over which
//    * properties of the forest to write. */
//   int write_treeid = 1;
//   int write_mpirank = 1;
//   int write_level = 1;
//   int write_element_id = 1;
//   int write_ghosts = 0;
//   t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, 0,
//                            0, num_data, &vtk_data);
//   T8_FREE (element_data);
// }

bool
fileExists (const std::string& filename)
{
  std::ifstream file (filename);
  return file.good ();
}

double
super_complex_function (double x, double y)
{
  // ------------------------------------------------------
  // 1. Fractal (Mandelbrot-like escape time simulation)
  double zx = 4.0 * x - 2.0;
  double zy = 4.0 * y - 2.0;
  double cx = zx, cy = zy;
  int max_iter = 100;
  int iter = 0;
  while (zx * zx + zy * zy < 4.0 && iter < max_iter) {
    double xtemp = zx * zx - zy * zy + cx;
    zy = 2.0 * zx * zy + cy;
    zx = xtemp;
    iter++;
  }
  double fractal = (double) iter / max_iter;

  // ------------------------------------------------------
  // 2. Sharp jump using sigmoid and singularity
  double jump = 0.0;
  if (fabs (1.0 - y * y) > 1e-6 && fabs (tan (y)) > 1e-6) {
    double exp_part = exp (3.0 * x * x);
    jump = sin (exp_part / tan (y)) * cos (1.0 / (1.0 - y * y)) * exp (y);
  }

  // ------------------------------------------------------
  // 3. Ripple pattern (radial sine decay)
  double dx = x - 0.5;
  double dy = y - 0.5;
  double r = sqrt (dx * dx + dy * dy);
  double theta = atan2 (dy, dx);
  double ripple = sin (25.0 * r - 5.0 * theta) * exp (-10.0 * r);

  // ------------------------------------------------------
  // 4. Swirl (vortex-like decay)
  double swirl = sin (10.0 * theta + 6.0 * r) * exp (-8.0 * r);

  // ------------------------------------------------------
  // 5. Turbulence (high-frequency oscillation)
  double turbulence = sin (60.0 * x) * cos (70.0 * y) * 0.2;

  // ------------------------------------------------------
  // Combine all components
  double result = 0.3 * fractal + 0.25 * jump + 0.2 * ripple + 0.15 * swirl + 0.2 * turbulence;

  return result;
}

float
calculate_correlation (float* array1, float* array2, int length)
{
  float sum1 = 0.0f, sum2 = 0.0f;
  float sum1_square = 0.0f, sum2_square = 0.0f;
  float product_sum = 0.0f;

  // Loop through the arrays and compute the necessary sums
  for (int i = 0; i < length; i++) {
    sum1 += array1[i];
    sum2 += array2[i];
    sum1_square += array1[i] * array1[i];
    sum2_square += array2[i] * array2[i];
    product_sum += array1[i] * array2[i];
  }

  // Calculate Pearson correlation coefficient
  float numerator = length * product_sum - sum1 * sum2;
  float denominator = sqrt ((length * sum1_square - sum1 * sum1) * (length * sum2_square - sum2 * sum2));

  if (denominator == 0) {
    return 0.0f;  // Avoid division by zero
  }

  return numerator / denominator;
}

// Function to calculate the indices that fall within the given bounds
void
find_indices (double xmin, double xmax, double ymin, double ymax, int** x_indices, int** y_indices, int* count)
{
  // Calculate the range of x indices that satisfy xmin <= xa[i] <= xmax
  int start_x = (int) ((xmin) / 0.3);
  int end_x = (int) ((xmax) / 0.3);

  // Calculate the range of y indices that satisfy ymin <= ya[i] <= ymax
  int start_y = (int) ((ymin + 90.0) / 0.3);
  int end_y = (int) ((ymax + 90.0) / 0.3);

  // Ensure that the indices are within valid bounds
  start_x = fmax (0, start_x);
  end_x = fmin (1201 - 1, end_x);
  start_y = fmax (0, start_y);
  end_y = fmin (601 - 1, end_y);

  // Allocate memory for the indices
  *count = (end_x - start_x + 1) * (end_y - start_y + 1);
  *x_indices = (int*) malloc (*count * sizeof (int));
  *y_indices = (int*) malloc (*count * sizeof (int));

  int idx = 0;
  // Collect all indices within the valid range
  for (int i = start_x; i <= end_x; ++i) {
    for (int j = start_y; j <= end_y; ++j) {
      (*x_indices)[idx] = i;
      (*y_indices)[idx] = j;
      ++idx;
    }
  }
}

// Function to create the combined index array and the corresponding coordinate array
void
create_combinations_and_coordinates (int* x_indices, int* y_indices, int count, double* xa, double* ya,
                                     int** index_combinations, double** coordinates)
{
  *index_combinations = (int*) malloc (count * 2 * sizeof (int));  // 2 * count (for x and y indices)
  *coordinates = (double*) malloc (count * 3 * sizeof (double));   // 3 * count (for x, y, z coordinates)

  int idx = 0;
  for (int i = 0; i < count; ++i) {
    (*index_combinations)[2 * idx] = x_indices[i];      // x_index
    (*index_combinations)[2 * idx + 1] = y_indices[i];  // y_index
    (*coordinates)[3 * idx] = xa[x_indices[i]];         // x value
    (*coordinates)[3 * idx + 1] = ya[y_indices[i]];     // y value
    (*coordinates)[3 * idx + 2] = 0.0;                  // z value (always 0)

    ++idx;
  }
}

// Function to iterate over the points where is_inside[i] == 1
void
print_matching_points (int* is_inside, double* coordinates, int num_points)
{
  printf ("Matching points (where is_inside[i] == 1):\n");

  // Iterate over all points
  for (int i = 0; i < num_points; ++i) {
    if (is_inside[i] == 1) {
      // Each point has 3 values: x, y, z (so we access coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2])
      printf ("[%.2f, %.2f, %.2f]\n", coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
    }
  }
}

void
fill_zvalues_arr_for_Correlation_wf_spline (t8_forest_t forest, float fileData_wf[721801], double* xa, double* ya)
{
  struct adapt_data_1d_wf_spline* adapt_data = (struct adapt_data_1d_wf_spline*) t8_forest_get_user_data (forest);
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme* scheme = t8_forest_get_scheme (forest);
  const t8_element_t* element;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    //t8_global_productionf ("test innen zwei \n");
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */
      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double volume = t8_forest_element_volume (forest, itree, element);
      double verts[3][3] = { 0 };
      t8_forest_element_coordinate (forest, itree, element, 0, verts[0]);
      t8_forest_element_coordinate (forest, itree, element, 1, verts[1]);
      t8_forest_element_coordinate (forest, itree, element, 2, verts[2]);
      double x_min, x_max, y_min, y_max;
      x_min = min (verts[0][0], verts[1][0]);
      x_min = min (x_min, verts[2][0]);
      x_max = max (verts[0][0], verts[1][0]);
      x_max = max (x_max, verts[2][0]);
      y_min = min (verts[0][1], verts[1][1]);
      y_min = min (y_min, verts[2][1]);
      y_max = max (verts[0][1], verts[1][1]);
      y_max = max (y_max, verts[2][1]);

      int *x_indices = NULL, *y_indices = NULL;
      int count = 0;
      const double tolerance = 1e-14;
      find_indices (x_min, x_max, y_min, y_max, &x_indices, &y_indices, &count);
      int* is_inside_arr;
      is_inside_arr = T8_ALLOC (int, count);
      // Arrays for combinations and coordinates
      int* index_combinations = NULL;
      double* coordinates = NULL;

      // Create combinations and coordinates
      create_combinations_and_coordinates (x_indices, y_indices, count, xa, ya, &index_combinations, &coordinates);

      t8_forest_element_points_inside (forest, itree, element, coordinates, count, is_inside_arr, tolerance);

      // Iterate over all points
      for (int i = 0; i < count; ++i) {
        if (is_inside_arr[i] == 1) {
          // Each point has 3 values: x, y, z (so we access coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2])
          fileData_wf[index_combinations[2 * i] * 601 + (600 - index_combinations[2 * i + 1])]
            = (float) AuswertungSinglescale_wf_spline (forest, scheme, coordinates[3 * i], coordinates[3 * i + 1],
                                                       itree, ielement, current_index);
          //printf("[%.2f, %.2f, %.2f]\n", coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
        }
      }
      free (x_indices);
      free (y_indices);
      free (index_combinations);
      free (coordinates);
      T8_FREE (is_inside_arr);
    }
  }
}

void
fill_zvalues_arr_for_Correlation_wb_spline (t8_forest_t forest, float fileData_wb[721801], double* xa, double* ya)
{
  struct adapt_data_1d_wb_spline* adapt_data = (struct adapt_data_1d_wb_spline*) t8_forest_get_user_data (forest);
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme* scheme = t8_forest_get_scheme (forest);
  const t8_element_t* element;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    //t8_global_productionf ("test innen zwei \n");
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */
      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double volume = t8_forest_element_volume (forest, itree, element);
      double verts[3][3] = { 0 };
      t8_forest_element_coordinate (forest, itree, element, 0, verts[0]);
      t8_forest_element_coordinate (forest, itree, element, 1, verts[1]);
      t8_forest_element_coordinate (forest, itree, element, 2, verts[2]);
      double x_min, x_max, y_min, y_max;
      x_min = min (verts[0][0], verts[1][0]);
      x_min = min (x_min, verts[2][0]);
      x_max = max (verts[0][0], verts[1][0]);
      x_max = max (x_max, verts[2][0]);
      y_min = min (verts[0][1], verts[1][1]);
      y_min = min (y_min, verts[2][1]);
      y_max = max (verts[0][1], verts[1][1]);
      y_max = max (y_max, verts[2][1]);

      int *x_indices = NULL, *y_indices = NULL;
      int count = 0;
      const double tolerance = 1e-14;
      find_indices (x_min, x_max, y_min, y_max, &x_indices, &y_indices, &count);
      int* is_inside_arr;
      is_inside_arr = T8_ALLOC (int, count);
      // Arrays for combinations and coordinates
      int* index_combinations = NULL;
      double* coordinates = NULL;

      // Create combinations and coordinates
      create_combinations_and_coordinates (x_indices, y_indices, count, xa, ya, &index_combinations, &coordinates);

      t8_forest_element_points_inside (forest, itree, element, coordinates, count, is_inside_arr, tolerance);

      // Iterate over all points
      for (int i = 0; i < count; ++i) {
        if (is_inside_arr[i] == 1) {
          // Each point has 3 values: x, y, z (so we access coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2])
          fileData_wb[index_combinations[2 * i] * 601 + (600 - index_combinations[2 * i + 1])]
            = (float) AuswertungSinglescale_wb_spline (forest, scheme, coordinates[3 * i], coordinates[3 * i + 1],
                                                       itree, ielement, current_index);
          //printf("[%.2f, %.2f, %.2f]\n", coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
        }
      }
      free (x_indices);
      free (y_indices);
      free (index_combinations);
      free (coordinates);
      T8_FREE (is_inside_arr);
    }
  }
}

// 5. Main program entry point
int
main (int argc, char** argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  cout << setprecision (15);
  cout << "File data is\n";
  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /* Creation of a basic two dimensional cmesh. */
  cout << "File data is\n";
  InitialisiereKoeff (p_mra, M0, M1, M2, M3, N0, N1, N2, N3);
  //Here you choose the correct coarse mesh.
  //t8_cmesh_t cmesh = t8_cmesh_new_basic (comm);
  t8_cmesh_t cmesh = t8_cmesh_new_debugging (comm);
  //t8_cmesh_t cmesh = t8_cmesh_new_octagon (comm);
  //t8_cmesh_t cmesh = t8_cmesh_new_complex_polygonal_shape (comm);
  //t8_cmesh_t cmesh = t8_cmesh_new_l_shape (comm);
  cout << "File data is\n";
  //t8_cmesh_t cmesh = t8_cmesh_new_earth (comm);
  t8_cmesh_vtk_write_file (cmesh, "cmesh_earth");
  const t8_scheme* scheme = t8_scheme_new_default ();
  //t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  //t8_forest_write_vtk (forest, "forest_lev0");
  cout << "File data is\n";
  int max_level = 8;
  levelgrid_map<t8_data_per_element_waveletfree_1d_gh>* grid
    = new levelgrid_map<t8_data_per_element_waveletfree_1d_gh> (max_level);
  levelgrid_map<t8_data_per_element_1d_gh>* grid_wb = new levelgrid_map<t8_data_per_element_1d_gh> (max_level);
  //forest = t8_forest_new_uniform (cmesh, scheme, 5, 0, comm);
  //t8_forest_write_vtk (forest, "forest_lev5");
  t8_forest_t forest;
  t8_forest_t forest_wb;
  t8_gloidx_t global_num_elements;
  //Write the data to a CSV file for later use in Python
  //This part writes the csv with the L1 and Linf error for the reference scheme.
  // Open the output CSV file
  //     std::ofstream outfile("error_ref_f3_octagon_p4.csv");
  //
  //     // Check if the file opened successfully
  //     if (!outfile.is_open()) {
  //         std::cerr << "Error opening file.\n";
  //         return 1;
  //     }
  //
  //   // Set precision for double values
  //   outfile << std::fixed << std::setprecision(15);
  //   outfile << "Level,Num Elements,L1 Error,Linf Error,p,C_thr\n";
  //   //outfile << "Level,Num Elements,L1 Error,L2 Error,Linf Error\n"; // CSV header
  // //test Plotting
  //   double C_thr=1.;
  //   int p=p_mra;
  // for(int level=1;level<=max_level;level++){
  //   forest=t8_create_init_mra_forest_wb_1D_func (grid_wb,comm,cmesh,scheme,1.,C_thr, 1.0,level,f3, 10,max_level);
  //   //calculate_rescale_wb_1D_func(forest);
  //   global_num_elements = t8_forest_get_global_num_elements (forest);
  //   const char* err_type;
  //   err_type = "L1";
  //   double l1_err=ErrorSinglescale(forest,10,err_type);
  //   err_type = "Linf";
  //   double linf_err=ErrorSinglescale(forest,10,err_type);
  //   //t8_global_productionf (" error: %f\n", err);
  //   // Write the data row to the CSV
  // outfile << level << ","
  //         << global_num_elements << ","
  //         << l1_err << ","
  //         << linf_err << ","
  //         << p << ","
  //         << C_thr << "\n";
  //   struct adapt_data_1d_wb_func *adapt_data;
  //   adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  //   sc_array_destroy (adapt_data->element_data);
  //   t8_global_productionf ("hier 2.\n");
  //   T8_FREE (adapt_data);
  //   //t8_forest_unref (&forest);
  // }
  //   outfile.close();

  //here is the part that writes the L1 and Linf error for the wavelet-based method
  double C_thr = 1.;
  int p = p_mra;
  t8_gloidx_t global_num_elements_next;
  // std::ofstream outfile("error_wb_f3_debug_p2.csv");
  // outfile << std::fixed << std::setprecision(15);
  //    outfile << "Level,Num Elements,L1 Error,Linf Error,p,C_thr\n";
  //    for(int level=1;level<=max_level;level++){
  //      grid_wb->erase_all();
  //      //C_thr=pow(2,-(level*(p_mra)));
  //      forest_wb=t8_create_init_mra_forest_wb_1D_func (grid_wb,comm,cmesh,scheme,1.,C_thr, 1.0,level,f3, 10,level);
  //      calculate_rescale_wf_1D_func(forest_wb);
  //      for(int i=0;i<level;i++){
  //        global_num_elements = t8_forest_get_global_num_elements (forest_wb);
  //        forest_wb=t8_thresholding_adapt(forest_wb);
  //        global_num_elements_next = t8_forest_get_global_num_elements (forest_wb);
  //        if(global_num_elements_next==global_num_elements){
  //          break;
  //        }
  //      }
  //      global_num_elements = t8_forest_get_global_num_elements (forest_wb);
  //      const char* err_type;
  //      err_type = "L1";
  //      double l1_err=ErrorSinglescale(forest_wb,10,err_type);
  //      err_type = "Linf";
  //      double linf_err=ErrorSinglescale(forest_wb,10,err_type);
  //      outfile << level << ","
  //            << global_num_elements << ","
  //            << l1_err << ","
  //            << linf_err << ","
  //            << p << ","
  //            << C_thr << "\n";
  //      struct adapt_data_1d_wb_func *adapt_data;
  //      adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_wb);
  //      sc_array_destroy (adapt_data->element_data);
  //      T8_FREE (adapt_data);
  //      //t8_forest_unref (&forest);
  //    }
  //      outfile.close();
  //
  //
  //
  //
  //
  //
  // //here is the part that writes the L1 and Linf error for the wavelet-free method
  // std::ofstream outfile2("error_wf_f7_octagon_p1.csv");
  // outfile2 << std::fixed << std::setprecision(15);
  //    outfile2 << "Level,Num Elements,L1 Error,Linf Error,p,C_thr\n";
  //       for(int level=1;level<=max_level;level++){
  //         grid->erase_all();
  //         //C_thr=pow(2,-(level*(p_mra)));
  //         forest=t8_create_init_mra_forest_wf_1D_func (grid,comm,cmesh,scheme,p_mra+1.,C_thr, 1.0,level,f7, 10,level);
  //         calculate_rescale_wf_1D_func(forest);
  //         for(int i=0;i<level;i++){
  //           global_num_elements = t8_forest_get_global_num_elements (forest);
  //           forest=t8_thresholding_adapt_wf(forest);
  //           global_num_elements_next = t8_forest_get_global_num_elements (forest);
  //            if(global_num_elements_next==global_num_elements){
  //              break;
  //            }
  //         }
  //         global_num_elements = t8_forest_get_global_num_elements (forest);
  //         const char* err_type;
  //         err_type = "L1";
  //         double l1_err=ErrorSinglescale_wf(forest,10,err_type);
  //         err_type = "Linf";
  //         double linf_err=ErrorSinglescale_wf(forest,10,err_type);
  //         outfile2 << level << ","
  //               << global_num_elements << ","
  //               << l1_err << ","
  //               << linf_err << ","
  //               << p << ","
  //               << C_thr << "\n";
  //         struct adapt_data_1d_wf_func *adapt_data;
  //         adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  //         sc_array_destroy (adapt_data->element_data);
  //         T8_FREE (adapt_data);
  //         //t8_forest_unref (&forest);
  //       }
  //         outfile2.close();

  // const size_t dataSize = 721801;
  // float fileData[dataSize];
  // std::string filename = "orig_data.h2o.lev.1000.bin";
  // if (fileExists(filename)) {
  //     std::cout << "File exists!" << std::endl;
  // } else {
  //     std::cerr << "File not found: " << filename << std::endl;
  // }
  // read_mptrac_data(filename, fileData, dataSize);

  //read the binary file for the MPTRAC data
  //read the binary file for the MPTRAC data
  //     cout << "davor?\n";
  //     float fileData[721801];
  //     ifstream inputFileStream("/home/veli/t8codeMA/t8code_build/tutorials/general/orig_data.h2o.lev.1097.bin", ios::in | ios::binary);
  //     //ifstream inputFileStream("/home/veli/t8codeMA/t8code_build/tutorials/general/orig_data.h2o.lev.1000-1747265754939.bin", ios::in | ios::binary);
  //     // "orig_data.h2o.lev.1000.bin"
  //     //"/home/veli/t8codeMA/t8code_build/tutorials/general/orig_data.h2o.lev.1000.bin"
  //     if (!inputFileStream) {
  //     std::cerr << "Error opening file!" << std::endl;
  //     return 1;  // or exit
  // }
  //     cout << "File data is\n";
  //     inputFileStream.read((char*) &fileData, 721801*sizeof(float));
  //
  //     /* If you want to output the data array
  //     for (int i=0; i<721801; i++)
  //         cout << fileData[i] << ", " << endl;
  //     cout << "\nFinished reading\n";
  //     */
  //     inputFileStream.close();
  //
  //     // for (int i = 0; i < 100; i++) { // Print the first 10 values
  //     //         std::cout << fileData[i] << std::endl;
  //     //     }
  //     const size_t dataSize = 721801;
  //     // Initialize min and max values using the first element
  //     float minValue = fileData[0];
  //     float maxValue = fileData[0];
  //
  //     // Iterate through the fileData array to find the min and max values
  //     for (size_t i = 1; i < dataSize; ++i) {
  //         if (fileData[i] < minValue) {
  //             minValue = fileData[i];  // Update min value
  //         }
  //         if (fileData[i] > maxValue) {
  //             maxValue = fileData[i];  // Update max value
  //         }
  //     }
  //
  //     // Printing the minimum and maximum values
  //     std::cout << "Minimum value in the data: " << minValue << std::endl;
  //     std::cout << "Maximum value in the data: " << maxValue << std::endl;
  //
  //     //2d interpolation using gsl library gsl_interp2d_bicubic gsl_interp2d_bilinear
  //     const gsl_interp2d_type *T = gsl_interp2d_bicubic;//or gsl_interp2d_bicubic
  //     //const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  //     const size_t N = 100;             /* number of points to interpolate */
  //     const size_t nx = 1201; /* x grid points */
  //     const size_t ny = 601; /* y grid points */
  //     double xa[nx]; /* define unit square */
  //     for(int i = 0; i < nx; i++) {
  //       xa[i] = 0.3*i;
  //     }
  //     double ya[ny];
  //     for(int i = 0; i < ny; i++) {
  //       ya[i] = -90.0+0.3*i;
  //     }
  //     double *za = T8_ALLOC (double, nx * ny);
  //     gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  //     gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  //     gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  //
  //     /* set z grid values */
  //     for (int iy = 600, i=0; iy>=0&&i<ny; --iy,++i){
  //       for (int ix = 0; ix < nx; ++ix){
  //         gsl_spline2d_set(spline, za, ix, i, fileData[ix*ny+iy]);
  //       }
  //     }
  //
  //     /* initialize interpolation */
  //     gsl_spline2d_init(spline, xa, ya, za, nx, ny);
  //     t8_global_productionf (" passt: \n");//+p_mra
  //     forest=t8_create_init_mra_forest_wf_1D_spline (grid,comm,cmesh,scheme,1.,0.0000001, 1.0,max_level,spline,xacc,yacc, 10,max_level);
  //     calculate_rescale_wf_1D_spline(forest);
  //     t8_mra_wf_write_vtu_spline ( forest,  "era5testref");
  //     float fileData_wf[721801];
  //     fill_zvalues_arr_for_Correlation_wf_spline(forest,fileData_wf,xa, ya);
  //     float correlation;
  //     correlation= calculate_correlation(fileData, fileData_wf, 721801);
  //     printf("Correlation: %f\n", correlation);
  //     for(int i=0;i<max_level;i++){
  //       //forest=t8_bottom_up_adapt(forest);
  //       forest=t8_thresholding_adapt_wf_spline(forest);
  //       //forest=t8_bottom_up_adapt(forest);
  //     }
  //     t8_mra_wf_write_vtu_spline ( forest,  "era5testwf");
  //
  //     fill_zvalues_arr_for_Correlation_wf_spline(forest,fileData_wf,xa, ya);
  //     correlation= calculate_correlation(fileData, fileData_wf, 721801);
  //     printf("Correlation after thresholding: %f\n", correlation);

  // forest_wb=t8_create_init_mra_forest_wb_1D_spline (grid_wb,comm,cmesh,scheme,1.+p_mra,0.01, 1.0,max_level,spline,xacc,yacc, 10,max_level);
  // calculate_rescale_wb_1D_spline(forest_wb);
  // float fileData_wb[721801];
  // //fill_zvalues_arr_for_Correlation_wb_spline(forest_wb,fileData_wb,xa, ya);
  // //correlation= calculate_correlation(fileData, fileData_wb, 721801);
  // printf("Correlation: %f\n", correlation);
  // for(int i=0;i<max_level;i++){
  //   //forest=t8_bottom_up_adapt(forest);
  //   forest_wb=t8_thresholding_adapt_wb_spline(forest_wb);
  //   //forest=t8_bottom_up_adapt(forest);
  // }
  // t8_mra_wb_write_vtu_spline ( forest_wb,  "era5testwb");
  //
  // //fill_zvalues_arr_for_Correlation_wb_spline(forest_wb,fileData_wb,xa, ya);
  // //correlation= calculate_correlation(fileData, fileData_wb, 721801);
  // printf("Correlation after thresholding: %f\n", correlation);
  //

  // // Define the bounds for the search
  // double xmin = 10.0, xmax = 50.0, ymin = -30.0, ymax = 30.0;
  //
  // // Containers to store the resulting indices
  // int *x_indices = NULL, *y_indices = NULL;
  // int count = 0;
  //
  // // Find indices within the given bounds
  // find_indices(xmin, xmax, ymin, ymax, xa, ya, nx, ny, &x_indices, &y_indices, &count);
  //
  // // Arrays for combinations and coordinates
  // int *index_combinations = NULL;
  // double *coordinates = NULL;
  //
  // // Create combinations and coordinates
  // create_combinations_and_coordinates(x_indices, y_indices, count, xa, ya,
  //                                     &index_combinations, &coordinates);
  //
  //                                     // Free dynamically allocated memory
  // free(x_indices);
  // free(y_indices);
  // free(index_combinations);
  // free(coordinates);

  forest = t8_create_init_mra_forest_wf_1D_func (grid, comm, cmesh, scheme, 1., 1., 1.0, max_level, f4, 10, max_level);
  calculate_rescale_wf_1D_func (forest);
  //const char* err_type_wb = "L2";
  //double err=ErrorSinglescale(forest,10,err_type);
  //err=ErrorSinglescale_wf(forest,10,err_type_wb);
  //t8_global_productionf (" error: %f\n", err);
  for (int i = 0; i < max_level; i++) {
    // forest = t8_bottom_up_adapt (forest);
    forest = t8_thresholding_adapt_wf (forest);
    // forest = t8_bottom_up_adapt (forest);
  }
  //  err=ErrorSinglescale_wf(forest,10,err_type);
  // t8_global_productionf (" error: %f\n", err);
  // // for(int i=0;i<max_level;i++){
  // //   //forest=t8_bottom_up_adapt(forest);
  // //   forest=t8_thresholding_adapt(forest);
  // // }
  t8_mra_wf_write_vtu (forest, "beforeprediction");
  set_min_level_wf (forest);
  add_old_level_element_data_wf (forest);
  int min_level = forest_get_min_level (forest);
  t8_global_productionf (" min_level: %i\n", min_level);
  for (int i = min_level; i < max_level - 1; i++) {
    //forest=t8_bottom_up_adapt(forest);
    t8_global_productionf (" i: %i\n", i);
    global_num_elements = t8_forest_get_global_num_leaf_elements (forest);
    forest = t8_prediction_adapt_wf (forest);
    global_num_elements_next = t8_forest_get_global_num_leaf_elements (forest);
    if (global_num_elements_next == global_num_elements) {
      break;
    }
  }
  t8_mra_wf_write_vtu (forest, "afterprediction");

  //
  //
  // forest_wb=t8_create_init_mra_forest_wb_1D_func (grid_wb,comm,cmesh,scheme,1.,0.1, 1.,max_level,f2, 10,max_level);
  // calculate_rescale_wb_1D_func(forest_wb);
  // //const char* err_type = "L2";
  // //double err=ErrorSinglescale(forest,10,err_type);
  // //double err=ErrorSinglescale(forest_wb,10,err_type);
  // //t8_global_productionf (" error: %f\n", err);
  // for(int i=0;i<max_level;i++){
  //   //forest=t8_bottom_up_adapt(forest);
  //   forest_wb=t8_thresholding_adapt(forest_wb);
  //   //forest_wb=t8_bottom_up_adapt(forest_wb);
  // }
  // //  err=ErrorSinglescale_wf(forest,10,err_type);
  // // t8_global_productionf (" error: %f\n", err);
  // // // for(int i=0;i<max_level;i++){
  // // //   //forest=t8_bottom_up_adapt(forest);
  // // //   forest=t8_thresholding_adapt(forest);
  // // // }
  // t8_mra_wb_write_vtu ( forest_wb,  "forestbeforegrading");
  // //t8_forest_write_vtk (forest, "prefix");
  // //add_old_level_element_data( forest);
  // forest_wb=t8_mra_balance ( forest_wb,1);
  //  t8_forest_write_vtk (forest_wb, "forestaftergrading");
  // const char *prefix_forest = "t8_mra_forest";
  // t8_forest_write_vtk (forest, prefix_forest);
  /* Write the adapted forest to a vtu file */

  //Run Prediction
  // forest=t8_create_init_mra_forest_wf_1D_func (grid,comm,cmesh,scheme,1.,1., 1.0,max_level,f3, 10,max_level);
  // calculate_rescale_wf_1D_func(forest);
  // for(int i=0;i<max_level;i++){
  //   forest=t8_thresholding_adapt_wf(forest);
  // }
  // add_old_level_element_data_wf(forest);

  struct adapt_data_1d_wb_func* adapt_data;
  adapt_data = (struct adapt_data_1d_wb_func*) t8_forest_get_user_data (forest);
  //
  // t8_global_productionf (" data test u%f.\n", t8_element_get_value (adapt_data, 0).u_coeff[0]);
  // t8_global_productionf (" data test d%f.\n", t8_element_get_value (adapt_data, 0).d_coeff[0]);
  // t8_global_productionf (" data test first %i.\n", t8_element_get_value (adapt_data, 0).first);
  // t8_global_productionf (" data test first %i.\n", t8_element_get_value (adapt_data, 0).second);
  // t8_global_productionf (" data test third %i.\n", t8_element_get_value (adapt_data, 0).third);
  // levelgrid_map<t8_data_per_element_1d_gh> grid(4);
  // //levelgrid_map<t8_data_per_element_1d_gh>* grid = new levelgrid_map<t8_data_per_element_1d_gh>(4);
  // struct t8_data_per_element_1d * element_data;
  // element_data=t8_create_element_data (grid, forest,F, 10, 4);
  // t8_forest_t adapt_forest=t8_bottom_up_init_adapt_forest (forest,0.00001,2.0,&grid,1);
  // const char *prefix_forest_adapt = "t8_mra_forest_adapt";
  // t8_forest_write_vtk (adapt_forest, prefix_forest_adapt);
  // 5.5. Return a value indicating successful execution
  /* Finalize the sc library */
  /* Destroy the forest. */
  /* Now you could continue working with the forest. */
  t8_global_productionf ("hier.\n");
  sc_array_destroy (adapt_data->element_data);
  t8_global_productionf ("hier 2.\n");
  T8_FREE (adapt_data);
  //T8_FREE ((adapt_data->element_data));
  //T8_FREE (adapt_data->element_data);
  t8_global_productionf ("vor unref.\n");
  grid->erase_all ();
  delete grid;
  t8_global_productionf ("vor unref refcount: %i\n", forest->rc.refcount);
  t8_forest_unref (&forest);

  ///////ab hier hatte ich auskommentiert

  //t8_cmesh_destroy (&cmesh);
  // Do something with the grid
  t8_global_productionf (" [step5] Destroyed forest.\n");
  // T8_FREE (element_data);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;  // Return 0 to indicate successful execution
}
