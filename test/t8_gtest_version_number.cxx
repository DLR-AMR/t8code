/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <gtest/gtest.h>
#include <t8_version.h>

/* Test whether the return value of t8_get_package_string
 * matches T8_PACKAGE_STRING.
 */
TEST (t8_gtest_version_number, getter_function)
{
  const char         *version_string = t8_get_package_string ();

  EXPECT_STREQ (version_string, T8_PACKAGE_STRING);
}

/* The version_string must have the format
 * "t8 version_number"
 */
TEST (t8_gtest_version_number, check_format_of_version_string)
{
  const char         *version_string = t8_get_package_string ();

  /* Check that version_string is not NULL.
   * This is an assertion since we cannot continue if it is NULL. */
  ASSERT_NE (version_string, nullptr);

  /* Copy it to a non-const char so that we can modify it by strtok */
  char                version_string_copy[BUFSIZ];
  strncpy (version_string_copy, version_string, BUFSIZ - 1);

  /* Check that the first three chars are 't8 ' */
  EXPECT_EQ (version_string_copy[0], 't');
  EXPECT_EQ (version_string_copy[1], '8');
  EXPECT_EQ (version_string_copy[2], ' ');

  /* Check that version_string == "t8 version_number" */
  const char         *version_number = t8_get_version_number ();
  EXPECT_STREQ (version_string + 3, version_number);
}

/* The version number should be a string "X.Y.Z-Something"
 * with X being the major and Y being the minor number and Z being the patch number. */
TEST (t8_gtest_version_number, check_version_number_has_major_minor_patch)
{
  const char         *version_number = t8_get_version_number ();
  const int           major = t8_get_version_major ();
  const int           minor = t8_get_version_minor ();
  const int           patch = t8_get_version_patch ();

  /* Copy version_number to non-const so that we can modify it by strtok */
  char                version_number_copy[BUFSIZ];
  strncpy (version_number_copy, version_number, BUFSIZ - 1);

  const char         *major_string = strtok (version_number_copy, ".");
  const char         *minor_string = strtok (NULL, ".");
  char               *point_string = strtok (NULL, ".");

  /* They should not be nullptr.
   * If they are, version_number does not contain two '.' */
  ASSERT_STRNE (major_string, nullptr);
  ASSERT_STRNE (minor_string, nullptr);
  ASSERT_STRNE (point_string, nullptr);

  /* Point string should match t8_get_version_point */
  const char         *point_string_check = t8_get_version_point_string ();
  EXPECT_STREQ (point_string, point_string_check);

  /* Extract patch number from point string.
   * The point string is either
   * "patch_number"
   * or
   * "patch_number-HASH"
   */
  const char         *patch_string = strtok (point_string, "-");

  /* Should not be nullptr */
  ASSERT_STRNE (patch_string, nullptr);

  /* Convert major, minor  and patch  from string to int */
  const int           major_string_to_int = atoi (major_string);
  const int           minor_string_to_int = atoi (minor_string);
  const int           patch_string_to_int = atoi (patch_string);

  /* Check that these match the t8code major and minor version number. */
  ASSERT_EQ (major_string_to_int, major);
  ASSERT_EQ (minor_string_to_int, minor);
  ASSERT_EQ (patch_string_to_int, patch);
}
