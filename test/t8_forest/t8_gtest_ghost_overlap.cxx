
#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_overlap.hxx>

/**
 * Unit Test for the ghost definition overlap
 */

TEST (t8_gtest_ghost_overlap, has_uniform_stretch_factor)
{
  std::array<double, 3> init_stretch_factor { 1.1, 1.4, 0.9 };
  t8_forest_ghost_definition_overlap ghost_definition (init_stretch_factor);

  EXPECT_TRUE (ghost_definition.has_uniform_stretch_factor ())
    << "The overlap ghost, init with a uniform stretch factor has no stretch factor.";

  std::array<double, 3> stretch_factor = ghost_definition.get_uniform_stretch_factor ();

  EXPECT_EQ (init_stretch_factor[0], stretch_factor[0])
    << "The uniform stretch factor given to the ghost definition differ at index 0 to the retrieved one.";
  EXPECT_EQ (init_stretch_factor[1], stretch_factor[1])
    << "The uniform stretch factor given to the ghost definition differ at index 1 to the retrieved one.";
  EXPECT_EQ (init_stretch_factor[2], stretch_factor[2])
    << "The uniform stretch factor given to the ghost definition differ at index 2 to the retrieved one.";
}

TEST (t8_gtest_ghost_overlap, set_uniform_stretch_factor)
{
  std::array<double, 3> stretch_factor { 1.1, 1.4, 0.9 };
  // t8_forest_ghost_definition_overlap* ghost_definition = new t8_forest_ghost_definition_overlap( );
  t8_forest_ghost_definition_overlap ghost_definition {};

  EXPECT_FALSE (ghost_definition.has_uniform_stretch_factor ())
    << "The overlap ghost, init with no uniform stretch factor has a uniform stretch factor.";

  ghost_definition.set_uniform_stretch_factor (stretch_factor);

  EXPECT_TRUE (ghost_definition.has_uniform_stretch_factor ())
    << "After set uniform stretch factor of ghost definition, no uniform stretch factor available.";

  // delete ghost_definition;
  // ghost_definition = nullptr;
}

TEST (t8_gtest_ghost_overlap, unable_uniform_stretch_factor)
{
  std::array<double, 3> stretch_factor { 1.1, 1.4, 0.9 };
  t8_forest_ghost_definition_overlap ghost_definition (stretch_factor);
  // t8_forest_ghost_definition_overlap* ghost_definition = new t8_forest_ghost_definition_overlap( stretch_factor );

  EXPECT_TRUE (ghost_definition.has_uniform_stretch_factor ())
    << "The overlap ghost, init with a uniform stretch factor has no stretch factor.";

  ghost_definition.unable_uniform_stretch_factor ();

  EXPECT_FALSE (ghost_definition.has_uniform_stretch_factor ())
    << "The overlap ghost, after unable the uniform stretch factor, has still a uniform stretch factor.";

  // delete ghost_definition;
  // ghost_definition = nullptr;
}
