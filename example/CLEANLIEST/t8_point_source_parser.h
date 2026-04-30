/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

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

// C++ standard library
#include<array>
#include<string>
#include<vector>

// t8code

/**
 * \brief This simple struct contains the properties of a H2 emission point source.
 */
struct PointSourceH2
{
  // Member variables
  std::string name {};                   ///< Some string identifier
  std::array<double, 3> coordinates {};  ///< The Cartesian coordinates
  double sourceRate {};                  ///< The rate at which H2 is "created" in the source
  double t_start {};                     ///< The point in time the source becomes active
  double t_end {};                       ///< The point in time the source becomes inactive
};

/**
 * \brief Read in point sources from a given CSV file and return them as std::vector.
 * 
 * \param [in] fileName 
 * \return std::vector containing all point sources read from \a fileName
 */
std::vector<PointSourceH2> t8_parse_point_sources_from_file(const std::string& fileName);
