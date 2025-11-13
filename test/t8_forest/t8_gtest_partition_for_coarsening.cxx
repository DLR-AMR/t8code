/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

// testing
#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

#include <t8_schemes/t8_default/t8_default.hxx>

// t8code
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_copy.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.h>


/** In this test we test the partition-for-coarsening functionality.
  * The test is based on the following main steps:
  *   - create and adapt an example forest and apply partition for coarsening
  *   - send the last 10 leaf elements of a processor to the next one
  *   - on each process, check whether a family is split across the border to the lower-rank process
  **/
// class t8_test_partition_for_coarsening_test: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>> {
class t8_test_partition_for_coarsening_test: public testing::TestWithParam<cmesh_example_base *> {

 public:
  
  /** Create an example forest, adapt it, and apply partition for coarsening. */
  // t8_forest_t
  void
  create_and_partition_example_forest ()
  {

    // For debugging purpose: Set to true to write vtk output.
    bool debug_write_vtk = true;

    // (1) Create cmesh and initially uniform forest.
    // ----------------------------------------------
    const int level = 3;
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;


    t8_cmesh_t mycmesh = tested_cmesh;


    // t8_cmesh_t other_cmesh;
    // t8_cmesh_init (&other_cmesh);

    // t8_cmesh_copy(other_cmesh,tested_cmesh, sc_MPI_COMM_SELF);

    // t8_cmesh_commit(other_cmesh, sc_MPI_COMM_SELF);

    // bool check = t8_cmesh_comm_is_valid(mycmesh,comm);
    // T8_ASSERT(check);

    // local_num_trees = t8_cmesh_get_num_local_trees (mycmesh);
    // global_num_trees = t8_cmesh_get_num_trees (mycmesh);

    t8_debugf("cmesh: Local number of trees: %i \n", t8_cmesh_get_num_local_trees (mycmesh));
    t8_debugf("cmesh: Global number of trees: %li \n", t8_cmesh_get_num_trees (mycmesh));

    // t8_productionf("Create uniform forest...");

    // // if(tested_cmesh->mpirank == 0) {
    int mpi_rank;
    sc_MPI_Comm_rank (comm, &mpi_rank);   

    // if(mpi_rank == 0) {
    // t8_forest_t forest = t8_forest_new_uniform (mycmesh, tested_scheme, level, 0, sc_MPI_COMM_WORLD);
    
    
    
    
    // //   tested_cmesh = t8_forest_get_cmesh(forest);

    // // mycmesh->

     t8_forest_t forest;
     t8_forest_init(&forest);

     t8_forest_set_cmesh(forest, mycmesh, comm);

     t8_forest_set_scheme(forest, tested_scheme);

     t8_forest_set_level(forest, level);

     t8_forest_commit(forest);



      t8_debugf("forest: Local number of leaves: %i \n", t8_forest_get_local_num_leaf_elements(forest));
      t8_debugf("forest: Global number of leaves: %li \n", t8_forest_get_global_num_leaf_elements(forest));
    // // // }
    // t8_productionf("... done!");

    // if (debug_write_vtk)
    // {
    //   std::string tmpName = "t8_test_pfc_uniform_";
    //   // tmpName = tmpName.append(t8_eclass_to_string[tested_eclass]);
    //   tmpName = tmpName.append(current_cmesh_name.c_str());
    //   t8_forest_write_vtk (forest, tmpName.c_str());
    // }

    // }

    // return forest;

  }

 protected:
  /** During SetUp, set the scheme and the eclass based on the current testing parameters.*/
  void
  SetUp () override
  {
    // const int scheme_id = std::get<0> (GetParam ());
    // tested_scheme = create_from_scheme_id (scheme_id);
    // tested_scheme->ref ();
    tested_scheme = t8_scheme_new_default();

    // cmesh_params::my_comms.pop_back();
    // cmesh_params::my_comms.push_back(sc_MPI_COMM_SELF);

    // size_t found = std::get<0> (GetParam ())->name.find (std::string ("empty"));
    // tested_cmesh = std::get<0> (GetParam ())->cmesh_create ();
    size_t found = (GetParam ())->name.find (std::string ("empty"));
    tested_cmesh = (GetParam ())->cmesh_create ();
    current_cmesh_name = (GetParam ())->name;
    
    if (found != std::string::npos) {
      /* Tests not working for empty cmeshes */
      t8_global_productionf("Skipping mesh: %s", current_cmesh_name.c_str() );
      GTEST_SKIP ();
    }
    found = GetParam ()->name.find (std::string ("bigmesh"));
    if (found != std::string::npos) {
      /* skip bigmeshes as they do not have vertices from which to build the connectivity */
      t8_global_productionf("Skipping mesh: %s", current_cmesh_name.c_str() );
      GTEST_SKIP ();
    }
    // found = GetParam ()->name.find (std::string ("hybrid"));
    // if (found != std::string::npos) {
    //   /* skip bigmeshes as they do not have vertices from which to build the connectivity */
    //   t8_global_productionf("Skipping mesh: %s", current_cmesh_name.c_str() );
    //   GTEST_SKIP ();
    // }

    if (t8_cmesh_get_dimension(tested_cmesh)<2) {
      t8_global_productionf("Skipping mesh: %s", current_cmesh_name.c_str() );
      GTEST_SKIP ();
    }

    t8_global_productionf("##############################################################################");
    t8_global_productionf("######## Testing mesh: %s", current_cmesh_name.c_str());
    t8_global_productionf("##############################################################################");

  }

  /** During TearDown, decrease the scheme's reference pointer to trigger its destruction.*/
  void
  TearDown () override
  {
    // t8_cmesh_unref (&tested_cmesh);
    // tested_scheme->unref ();
    t8_global_productionf("TearDown!");
  }

  // Member variables: The currently tested scheme and eclass.
  const t8_scheme *tested_scheme;
  // t8_eclass_t tested_eclass;
  t8_cmesh_t tested_cmesh;
  std::string current_cmesh_name;
};

// The test's main function.
TEST_P (t8_test_partition_for_coarsening_test, test_partition_for_coarsening)
{

  sc_MPI_Barrier (sc_MPI_COMM_WORLD);


  // Create an example forest and repartition it using partition for coarsening.
  // t8_forest_t forest = 
  create_and_partition_example_forest ();

  sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  sleep(1);


  // Memory deallocation: Destroy the forest, but keep the scheme.
  // tested_scheme->ref ();
  // t8_forest_unref (&forest);
}

// Instantiate parameterized test to be run for all schemes.
INSTANTIATE_TEST_SUITE_P (t8_gtest_partititon_for_coarsening, t8_test_partition_for_coarsening_test, 
                        AllCmeshsParam,
                          pretty_print_base_example);
                          
// Instantiate parameterized test to be run for all schemes.
// INSTANTIATE_TEST_SUITE_P (t8_gtest_partititon_for_coarsening, t8_test_partition_for_coarsening_test, 
//                         testing::Combine (AllSchemeCollections, AllCmeshsParam),
//                           pretty_print_base_example_scheme);

