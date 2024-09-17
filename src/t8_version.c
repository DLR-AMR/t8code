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

#include <t8_version.h>

const char *
t8_get_package_string ()
{
  return T8_PACKAGE_STRING;
}

const char *
t8_get_version_point_string ()
{
  return T8_VERSION_POINT_STRING;
}

const char *
t8_get_version_number ()
{
  return T8_VERSION;
}

int
t8_get_version_major ()
{
  return T8_VERSION_MAJOR;
}

int
t8_get_version_minor ()
{
  return T8_VERSION_MINOR;
}

int
t8_get_version_patch ()
{
#ifndef T8_CMAKE_BUILD
  const char *version_point = t8_get_version_point_string ();

  if (version_point == NULL || strlen (version_point) == 0) {
    /* The patch version is invalid. */
    t8_global_errorf ("ERROR: Point version string is NULL or empty.\n");
    /* We abort in debugging mode and continue with value -1 otherwise. */
    T8_ASSERT (version_point != NULL && strlen (version_point) > 1);
    return -1;
  }
  /* We expect the string to look like "X-HASH" or "X", 
   * with X an integer and HASH a string. */
  char *error_check_string;
  const long patch_number_long = strtol (version_point, &error_check_string, 10);
  /* strotl returns to error_check_string the first invalid character of versiont_point,
   * thus the first character that is not a number. */

  /* The first invalid character must not be the first character of version_point */
  if (*error_check_string == version_point[0]) {
    /* The string does not start with a number */
    t8_global_errorf ("ERROR: Point version string '%s' does not begin with patch number.\n", version_point);
    /* Abort in debugging mode, continue with -1 else. */
    T8_ASSERT (*error_check_string == version_point[0]);
    return -1;
  }

  /* Convert to integer */
  const int patch_number = patch_number_long;
  /* Double check conversion */
  T8_ASSERT (patch_number == patch_number_long);
  if (patch_number < 0) {
    t8_global_errorf ("ERROR: Patch number %i is not >=0\n", patch_number);
    /* Abort in debugging mode, continue with negative number else. */
    T8_ASSERT (patch_number >= 0);
  }

  /* Return the patch number */
  return patch_number;
#else
  return T8_VERSION_PATCH;
#endif
}
