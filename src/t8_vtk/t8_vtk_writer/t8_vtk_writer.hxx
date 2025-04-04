/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#ifndef T8_VTK_WRITER_HXX
#define T8_VTK_WRITER_HXX

#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include "t8_forest/t8_forest_types.h"
#include "t8_vtk/t8_vtk_writer_helper.hxx"
#include "t8_vtk/t8_vtk_write_ASCII.hxx"

#include <string>
#include <t8_vtk.h>
#include <t8_types/t8_vec.hxx>

/**
 * A class that controls the writing of vtk files for cmeshes or forests. 
 * 
 * \tparam grid_t can be a forest or a cmesh. 
 */
template <typename grid_t>
class vtk_writer {
 public:
  /**
   * Construct a new vtk writer object. All parameters are set to false by default. By default no data is used and
   * \a num_data is set to zero. A default \a fileprefix is NOT given. 
   * 
   * \param write_treeid True, if we want to write the tree id of every element.
   * \param write_mpirank True, if we want to write the mpirankof every element.
   * \param write_level True, if we want to write the level of every element. Uses level 0 if used for a cmesh.
   * \param write_element_id True, if we want to write the element id of every element. Ignored if used for a cmesh.
   * \param write_ghosts True, if we want to write the ghost elements, too. 
   * \param curved_flag True, if we want to use quadratic vtk cells. Uses the geometry of the grid to evaluate points between corners. 
   * \param fileprefix The prefix of the output-file.
   * \param num_data The number of data-fields to print.
   * \param data The data to use.
   * \param comm The communicator for parallel output.
   */
  vtk_writer (const bool write_treeid, const bool write_mpirank, const bool write_level, const bool write_element_id,
              const bool write_ghosts, const bool curved_flag, std::string fileprefix, const int num_data,
              t8_vtk_data_field_t *data, sc_MPI_Comm comm)
    : write_treeid (write_treeid), write_mpirank (write_mpirank), write_level (write_level),
      write_element_id (write_element_id), write_ghosts (write_ghosts), curved_flag (curved_flag),
      fileprefix (fileprefix), num_data (num_data), data (data), comm (comm)
  {
  }

  /**
   * Construct a new vtk writer object. All parameters are set to false.
   * 
   * \param[in] fileprefix 
   * \param[in] comm 
   */
  vtk_writer (std::string fileprefix, sc_MPI_Comm comm): fileprefix (fileprefix), comm (comm)
  {
  }

  /**
   * A vtk-writer function that uses the vtk API.
   * 
   * \param[in] grid The forest or cmesh that is translated.
   * \return true, if writing was successful. 
   * \return false if writing was not successful. 
   */
  bool
  write_with_API (const grid_t grid)
  {
    return write_vtk (grid);
  }

  /**
   * A vtk-writer function that uses the vtk API
   * 
   * \param[in] grid The forest or cmesh that is translated
   * \return true 
   * \return false 
   */
  bool
  write_ASCII (const grid_t grid);

  /**
   * Set the write treeid flag. Set to true, if you want to write the tree id of every element.
   * 
   * \param[in] write_treeid true or false
   */
  inline void
  set_write_treeid (const bool write_treeid)
  {
    this->write_treeid = write_treeid;
  }

  /**
   * Set the write mpirank flag. Set to true, if you want to write the mpirank of every element.
   * 
   * \param[in] write_mpirank true or false
   */
  inline void
  set_write_mpirank (const bool write_mpirank)
  {
    this->write_mpirank = write_mpirank;
  }

  /**
   * Set the write level flag. Set to true, if you want to write the level of every element.
   * 
   * \param[in] write_level true or false
   */
  inline void
  set_write_level (const bool write_level)
  {
    this->write_level = write_level;
  }

  /**
   * Set the write element id flag. Set to true, if you want to write the element id of every element.
   * 
   * \param[in] write_element_id true or false
   */
  inline void
  set_write_element_id (const bool write_element_id)
  {
    this->write_element_id = write_element_id;
  }

  /**
   * Set the write ghosts flag. Set to true, if you want to write the ghost elements, too.
   * 
   * \param[in] write_ghosts true or false
   */
  inline void
  set_write_ghosts (const bool write_ghosts)
  {
    this->write_ghosts = write_ghosts;
  }

  /**
   * Set the curved flag. Set to true, if you want to use quadratic vtk cells. 
   * Uses the geometry of the grid to evaluate points between corners.
   * 
   * \param[in] curved_flag true or false
   */
  inline void
  set_curved_flag (const bool curved_flag)
  {
    this->curved_flag = curved_flag;
  }

  /**
   * Set the fileprefix for the output files.
   * \param[in] fileprefix 
   */
  inline void
  set_fileprefix (std::string fileprefix)
  {
    this->fileprefix = fileprefix;
  }

 private:

  bool write_treeid = false;
  bool write_mpirank = false;
  bool write_level = false;
  bool write_element_id = false;
  bool write_ghosts = false;
  bool curved_flag = false;

protected:

  std::string fileprefix;
  int num_data = 0;
  t8_vtk_data_field_t *data = NULL;
  sc_MPI_Comm comm;

};

#endif /* T8_VTK_WRITER_HXX */
