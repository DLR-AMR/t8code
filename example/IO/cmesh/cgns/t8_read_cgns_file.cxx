/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/** \file t8_read_cgns_file.cxx
 * This example reads a CGNS file into a cmesh, prints the zones and boundary conditions that were
 * found and optionally writes the cmesh to a vtk file.
 *
 * Usage: t8_read_cgns_file -f mesh.cgns [-o output_prefix]
 */

#include <sc_options.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_cgns_reader.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>
#include <t8_vtk/t8_vtk_writer.h>

#include <cstring>

/** Read the file, print its content and write it to vtk.
 * \return 0 on success, 1 if the file could not be read.
 */
static int
t8_read_cgns_file (const char *filename, const char *vtk_prefix)
{
  t8_cmesh_cgns_reader reader (filename);
  t8_cmesh_t cmesh = reader.read (sc_MPI_COMM_WORLD);
  if (cmesh == nullptr) {
    t8_global_errorf ("Could not read %s.\n", filename);
    return 1;
  }

  t8_global_productionf ("Read a cmesh of dimension %d with %lld trees from %s.\n", reader.get_dimension (),
                         static_cast<long long> (t8_cmesh_get_num_trees (cmesh)), filename);
  t8_global_productionf ("Zones:\n");
  for (const t8_cmesh_cgns_zone &zone : reader.get_zones ()) {
    t8_global_productionf (
      "  %-32s trees %lld to %lld, %lld vertices\n", zone.name.c_str (), static_cast<long long> (zone.first_tree),
      static_cast<long long> (zone.first_tree + zone.num_trees - 1), static_cast<long long> (zone.num_vertices));
  }
  t8_global_productionf ("Boundary conditions:\n");
  for (const t8_cmesh_cgns_boundary_condition &bc : reader.get_boundary_conditions ()) {
    t8_global_productionf ("  %-32s type %s, family \"%s\", zone %s, %lld tree faces\n", bc.name.c_str (),
                           bc.type_name.c_str (), bc.family.c_str (), reader.get_zones ()[bc.zone].name.c_str (),
                           static_cast<long long> (bc.num_faces));
  }

  /* Count the tree faces at the domain boundary. */
  long long num_boundary_faces = 0;
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    for (int face = 0; face < t8_eclass_num_faces[eclass]; ++face) {
      if (t8_cmesh_tree_face_is_boundary (cmesh, itree, face)) {
        ++num_boundary_faces;
      }
    }
  }
  t8_global_productionf ("%lld tree faces lie on the domain boundary.\n", num_boundary_faces);

  if (strcmp (vtk_prefix, "") != 0) {
    if (t8_cmesh_vtk_write_file (cmesh, vtk_prefix)) {
      t8_global_productionf ("Wrote the cmesh to %s.pvtu.\n", vtk_prefix);
    }
    else {
      t8_global_errorf ("Could not write the cmesh to %s.\n", vtk_prefix);
    }
  }
  t8_cmesh_destroy (&cmesh);
  return 0;
}

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  int helpme = 0;
  const char *filename = nullptr;
  const char *vtk_prefix = nullptr;
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'f', "file", &filename, "", "The CGNS file to read.");
  sc_options_add_string (opt, 'o', "output", &vtk_prefix, "",
                         "If given, the cmesh is written to vtk files with this prefix.");
  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  int retval = 0;
  if (helpme) {
    t8_global_productionf ("Read a CGNS file into a cmesh and print its zones and boundary conditions.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed < 0 || strcmp (filename, "") == 0) {
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    retval = 1;
  }
  else {
    retval = t8_read_cgns_file (filename, vtk_prefix);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return retval;
}
