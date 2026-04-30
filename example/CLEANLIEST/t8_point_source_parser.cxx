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
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>

// t8code
#include "t8_point_source_parser.h"
#include <t8.h>


// Trim whitespace from start and end of a string
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    size_t end = str.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// Split a CSV line into fields, handling quoted fields
static std::vector<std::string> splitCSVLine(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quotes = false;
    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if (c == ',' && !in_quotes) {
            fields.push_back(trim(field));
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(trim(field)); // Add the last field
    return fields;
}



std::vector<PointSourceH2> t8_parse_point_sources_from_file(const std::string& fileName)
{
  
  // Define point-source vector
  std::vector<PointSourceH2> pointSourceVec;

  // Throw error if file can not be read.
  std::ifstream file(fileName);
  if (!file.is_open()) {
    SC_ABORTF ("Could not open the file %s.\n", fileName.c_str ());
  }

  // Prepare parsing.
  std::string line;
  // Skip header line
  std::getline(file, line);

  // Loop over lines
  while (std::getline(file, line)) {
      auto fields = splitCSVLine(line);
      if (fields.size() != 7) {
          std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
          continue;
      }
      try {
          // Read row data into point source.
          PointSourceH2 curSource;
          curSource.coordinates[0] = std::stod(fields[0]);
          curSource.coordinates[1] = std::stod(fields[1]);
          curSource.coordinates[2] = std::stod(fields[2]);
          curSource.name = fields[3];
          curSource.sourceRate = std::stod(fields[4]);
          curSource.t_start = std::stod(fields[5]);
          curSource.t_end = std::stod(fields[6]);

          // Add to point-source vector.
          pointSourceVec.push_back(curSource);

      } catch (const std::exception& e) {
          std::cerr << "Warning: Could not parse line: " << line << " (" << e.what() << ")" << std::endl;
      }
  }

  // Return point-source vector
  return pointSourceVec;
}