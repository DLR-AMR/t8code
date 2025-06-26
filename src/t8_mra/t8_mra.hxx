// t8 mra hxx
// ich brauche folgendes:
// element data hier ohne Details
// die ganzen FUnktionen wie vorher
// Auswertung MS
// AUswertung SS
//
//
// Bottom Up Init:
// input: cmesh, ctresh, gamma, grid hierarchy, soll forest direkt adapten
// und grid hierarchy Daten eintragen
// thresholding:
// Ctresh, grid hierarchy: Soll Gitterhierarchie und forest anpassen
// Gitterhierarchie nur signifikant setzen
// Prediction:
// grid hierarchy anpassen und forest adaptieren
// nehmen als Input immer forest und grid hierarchy
//
// Element data ohne Detail coefficients
//
// lmi Funktionen

#include <t8.h>                          /* General t8code header, always include this. */
#include <t8_cmesh.hxx>                  /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h> /* forest definition and basic interface. */
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_io.h>                                       /* forest io interface. */
#include <t8_schemes/t8_default/t8_default.hxx>                           /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>                              /* geometrical information */
#include "vecmat.hxx"
#include "basis_functions.hxx"
#include "mask_coefficients.hxx"
#include "dunavant.hxx"
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
#include "t8_unordered_dense.hxx"
#include "t8_gsl.hxx"
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <fstream>
#include <time.h>

mat M0, M1, M2, M3, N0, N1, N2, N3;
typedef double (*func) (double, double);
/* declaring the matrices */
/* defines the maximum refinement level (should currently be less than approx. 10 depending on how many basecells there are) */
typedef int8_t t8_dtri_cube_id_t;
// Precompute powers of 4 up to level 29(T8_DTRI_MAXLEVEL)
#define MAX_LEVEL T8_DTRI_MAXLEVEL
//double gamma_scaling[2][MAX_LEVEL + 1];
/* To store the data for the 3d case (e.g. in AuswertungSinglescale) */
struct double_3d_array
{
  double dim_val[3];
};

// /* This is our own defined data that we will pass on to the
//  * adaptation callback. */
// struct t8_step7_adapt_data
// {
//   t8_3D_point midpoint;             /* The midpoint of our sphere. */
//   double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
//   double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
//   sc_array_t *element_data;
// };

/* We can drop the volume of the element and the level, hasFather, haschilds, childs ids, father ids*/
/* This struct stores the element data */
struct t8_data_per_element_1d_predict
{
  double u_coeff[M_mra];  //single-scale coefficients for all dof/ basis polynomials
  uint64_t lmi;
  unsigned int level_before_predict : 5;
};

/* We can drop the volume of the element and the level, hasFather, haschilds, childs ids, father ids*/
/* This struct stores the element data */
struct t8_data_per_element_1d
{
  double u_coeff[M_mra];  //single-scale coefficients for all dof/ basis polynomials
  uint64_t lmi;
};

/* This struct stores the element data */
struct t8_data_per_element_3d
{
  double u_coeff_d1[M_mra];  //single-scale coefficients for all dof/ basis polynomials
  double u_coeff_d2[M_mra];
  double u_coeff_d3[M_mra];
  uint64_t lmi;
};

/* This struct stores the element data */
struct t8_data_per_element_3d_predict
{
  double u_coeff_d1[M_mra];  //single-scale coefficients for all dof/ basis polynomials
  double u_coeff_d2[M_mra];
  double u_coeff_d3[M_mra];
  uint64_t lmi;
  unsigned int level_before_predict : 5;
};

/* This struct stores the element data */
struct adapt_data_1d_wb_func
{
  double gamma;
  double C_thr;
  double c_rescale;
  int max_level;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_1d_gh> *grid_map_ptr;
  func my_func;
};

/* This struct stores the element data */
struct adapt_data_1d_wb_spline
{
  double gamma;
  double C_thr;
  double c_rescale;
  int max_level;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_1d_gh> *grid_map_ptr;
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
};

/* This struct stores the element data */
struct adapt_data_1d_wf_func
{
  double gamma;
  double C_thr;
  double c_rescale;
  int max_level;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_map_ptr;
  func my_func;
};

/* This struct stores the element data */
struct adapt_data_1d_wf_spline
{
  double gamma;
  double C_thr;
  double c_rescale;
  int max_level;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_map_ptr;
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
};

/* This struct stores the element data */
struct adapt_data_3d_wb_func
{
  double gamma;
  double C_thr;
  double c_rescale1;
  double c_rescale2;
  double c_rescale3;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_3d_gh> *grid_map_ptr;
  func my_func1;
  func my_func2;
  func my_func3;
};

/* This struct stores the element data */
struct adapt_data_3d_wf_func
{
  double gamma;
  double C_thr;
  double c_rescale1;
  double c_rescale2;
  double c_rescale3;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_waveletfree_3d_gh> *grid_map_ptr;
  func my_func1;
  func my_func2;
  func my_func3;
};

/* This struct stores the element data */
struct adapt_data_3d_wb_spline
{
  double gamma;
  double C_thr;
  double c_rescale1;
  double c_rescale2;
  double c_rescale3;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_3d_gh> *grid_map_ptr;
  // spline my_func1;
  // spline my_func2;
  // spline my_func3;
  func my_func1;
  func my_func2;
  func my_func3;
};
/* This struct stores the element data */
struct adapt_data_3d_wf_spline
{
  double gamma;
  double C_thr;
  double c_rescale1;
  double c_rescale2;
  double c_rescale3;
  int current_level;
  sc_array_t *element_data;
  levelgrid_map<t8_data_per_element_waveletfree_3d_gh> *grid_map_ptr;
  // spline my_func1;
  // spline my_func2;
  // spline my_func3;
  func my_func1;
  func my_func2;
  func my_func3;
};

struct t8_data_per_element_1d *
t8_create_element_data (levelgrid_map<t8_data_per_element_1d_gh> &grid_hierarchy, t8_forest_t forest, func F,
                        const int rule, const int max_lev);
/* Set the value of an element to a given entry. */
// void
//  t8_element_set_value (const adapt_data_1d_wo_gh *adapt_data, t8_locidx_t ielement, t8_data_per_element_1d_levelwise_bottom_up elem_data);
//  /* Set the value of an element to a given element data. */
// void
//  t8_element_set_element (const adapt_data_1d_wo_gh *adapt_data, t8_locidx_t ielement, t8_data_per_element_1d_levelwise_bottom_up element);
//  /* Get the value of an element. */
// t8_data_per_element_1d_levelwise_bottom_up
//  t8_element_get_value (const adapt_data_1d_wo_gh *adapt_data, t8_locidx_t ielement);
t8_forest_t
t8_create_init_mra_forest_wb_1D_func (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, sc_MPI_Comm comm,
                                      t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma, double C_thr,
                                      double c_rescale, int current_level, func my_func, const int rule, int max_level);
void
calculate_rescale_wb_1D_func (t8_forest_t forest);
void
t8_mra_wb_write_vtu (t8_forest_t forest, const char *prefix);
void
t8_mra_wf_write_vtu (t8_forest_t forest, const char *prefix);
void
calculate_rescale_wf_1D_func (t8_forest_t forest);
t8_forest_t
t8_create_init_mra_forest_wf_1D_spline (levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_hierarchy,
                                        sc_MPI_Comm comm, t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma,
                                        double C_thr, double c_rescale, int current_level, gsl_spline2d *spline,
                                        gsl_interp_accel *xacc, gsl_interp_accel *yacc, const int rule, int max_level);
t8_forest_t
t8_create_init_mra_forest_wb_1D_spline (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, sc_MPI_Comm comm,
                                        t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma, double C_thr,
                                        double c_rescale, int current_level, gsl_spline2d *spline,
                                        gsl_interp_accel *xacc, gsl_interp_accel *yacc, const int rule, int max_level);
void
add_old_level_element_data (t8_forest_t forest);
void
add_old_level_element_data_wf (t8_forest_t forest);
double
ErrorSinglescale (t8_forest_t forest, int rule, const char *err_type);
double
ErrorSinglescale_wf (t8_forest_t forest, int rule, const char *err_type);
t8_forest_t
t8_thresholding_adapt_wb_spline (t8_forest_t forest);
t8_forest_t
t8_bottom_up_adapt (t8_forest_t forest);
void
t8_mra_wf_write_vtu_spline (t8_forest_t forest, const char *prefix);
void
t8_mra_wb_write_vtu_spline (t8_forest_t forest, const char *prefix);
t8_forest_t
t8_thresholding_adapt_wf_spline (t8_forest_t forest);
t8_forest_t
t8_thresholding_adapt_spline (t8_forest_t forest);
t8_forest_t
t8_prediction_adapt_wf (t8_forest_t forest);
void
calculate_rescale_wb_1D_spline (t8_forest_t forest);
int
forest_get_min_level (t8_forest_t forest);
void
calculate_rescale_wf_1D_spline (t8_forest_t forest);
void
set_min_level_wf (t8_forest_t forest);
double
AuswertungSinglescale_wf (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                          t8_locidx_t ielement, t8_locidx_t current_index);
double
AuswertungSinglescale (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                       t8_locidx_t ielement, t8_locidx_t current_index);
double
AuswertungSinglescale_wb_spline (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                                 t8_locidx_t ielement, t8_locidx_t current_index);
double
AuswertungSinglescale_wf_spline (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                                 t8_locidx_t ielement, t8_locidx_t current_index);
t8_forest_t
t8_thresholding_adapt (t8_forest_t forest);
t8_forest_t
t8_mra_balance (t8_forest_t forest, int partition);
t8_forest_t
t8_thresholding_adapt_wf (t8_forest_t forest);
t8_forest_t
t8_create_init_mra_forest_wf_1D_func (levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_hierarchy,
                                      sc_MPI_Comm comm, t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma,
                                      double C_thr, double c_rescale, int current_level, func my_func, const int rule,
                                      int max_level);
t8_forest_t
t8_bottom_up_init_adapt_forest (t8_forest_t forest, double C_thr, double gamma,
                                levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, int max_lev);
