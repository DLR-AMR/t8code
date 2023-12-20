/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <t8.h>                                     /* General t8code header, always include this. */
#include <sc_options.h>                             /* CLI parser */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_cmesh_vtk_writer.h>
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                 /* save forest */
#include <t8_forest/t8_forest_geometrical.h>        /* geometrical information of the forest */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* Linear geometry calculation of trees */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>    /* Curved geometry calculation of trees */
#include <t8_cmesh_readmshfile.h>                                         /* msh file reader */
#include <string>                                                         /* std::string */
#include <array>                                                          /* std::array */
#include <unordered_map>
#include <utility>
#include <vector>

#define NUM_SPARKLES 60
#define NUM_CONSTANT_SPARKLES 30
#define LEVEL_SPARKLE 4
#define NUM_ROUNDS 50


  int first_changing_sparkle = 0;

const std::unordered_map<t8_gloidx_t, int> tree_map
  = { { 0, 2 },  { 1, 0 },  { 2, 0 },  { 3, 0 },  { 4, 2 },  { 5, 0 },  { 6, 3 },
      { 7, 2 },  { 8, 0 },  { 9, 3 },  { 10, 3 }, { 11, 0 }, { 12, 1 }, { 13, 2 },
      { 14, 0 }, { 15, 0 }, { 16, 1 }, { 17, 1 }, { 18, 0 }, { 19, 0 }, { 20, 0 } };

std::vector<std::pair<double, double>> sparkles;

int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                   t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int elem_level = ts->t8_element_level (elements[0]);
  const t8_gloidx_t gtreeid = t8_cmesh_get_global_id (t8_forest_get_cmesh (forest_from), which_tree);
  if (elem_level < tree_map.at (gtreeid)) {
    return 1;
  }
  else if (elem_level > tree_map.at (gtreeid) && is_family) {
    return -1;
  }
  return 0;
}

int
t8_adapt_sparcle_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_eclass_scheme_c *ts, const int is_family, const int num_elements,
                           t8_element_t *elements[])
{
  const int elem_level = ts->t8_element_level (elements[0]);
  const t8_gloidx_t gtreeid = t8_cmesh_get_global_id (t8_forest_get_cmesh (forest_from), which_tree);
  
  if (gtreeid == 18 || gtreeid == 19 || gtreeid == 17 || gtreeid == 16) {
    return 0;
  }
  
  int is_inside = 0;
  for (auto sparkle : sparkles) {
    const double sparkle_array[3] = { sparkle.first, sparkle.second, 0 };
    if (t8_forest_element_point_inside (forest_from, which_tree, *elements, sparkle_array, 0.01)) {
      is_inside = 1;
      break;
    }
  }
  if (elem_level < (is_inside ? LEVEL_SPARKLE + 1 : tree_map.at (gtreeid))) {
    return 1;
  }
  return 0;
}

int
t8_generate_sparkles ()
{
  sparkles.resize (first_changing_sparkle);
  for (int i_sparkle = first_changing_sparkle; i_sparkle < NUM_SPARKLES; ++i_sparkle) {
    std::pair<double, double> sparkle;
    sparkle.first = rand () / (double) RAND_MAX * 730.0 - 730 / 2.0;
    sparkle.second = rand () / (double) RAND_MAX * 730.0;
    sparkles.push_back (sparkle);
  }
  first_changing_sparkle = NUM_CONSTANT_SPARKLES;
  return 0;
}

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);
  t8_forest_t forest, forest_new;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file ("t8_christmas_card", 0, comm, 2, 0, 1);
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0, comm);
  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest, t8_adapt_callback, 1);
  t8_forest_commit (forest_new);
  forest = forest_new;
  t8_forest_write_vtk_ext (forest, "christmas_card_0", 1, 1, 1, 1, 0, 1, 0, 0, NULL);
  t8_cmesh_vtk_write_file (cmesh, "christmas_card_cmesh", 1.0);
  for (int rounds = 1; rounds < NUM_ROUNDS + 1; ++rounds) {
    t8_generate_sparkles ();
    t8_forest_init (&forest_new);
    t8_forest_set_adapt (forest_new, forest, t8_adapt_sparcle_callback, 1);
    //t8_forest_set_balance(forest_new, forest, 1);
    t8_forest_commit (forest_new);
    forest = forest_new;

    t8_forest_init (&forest_new);
    t8_forest_set_adapt (forest_new, forest, t8_adapt_callback, 0);
    //t8_forest_set_balance(forest_new, forest, 1);
    t8_forest_commit (forest_new);
    forest = forest_new;
    std::string filename = "christmas_card_" + std::to_string (rounds);
    t8_forest_write_vtk_ext (forest, filename.c_str (), 1, 1, 1, 1, 0, 1, 0, 0, NULL);
  }

  t8_forest_unref (&forest);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}