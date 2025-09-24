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
  return T8_VERSION_PATCH;
}
