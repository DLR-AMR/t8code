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

/**
 * \file Implementation for a gmsh-reader.
 * 
 */

#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <src/t8_IO/t8_reader/t8_gmsh_reader/t8_gmsh_reader.hxx>
#include <src/t8_cmesh/t8_cmesh_types.h>
#include <src/t8_cmesh/t8_cmesh_stash.h>

/* The supported number of gmesh tree classes.
 * Currently, we only support first order trees.
 */
#define       T8_NUM_GMSH_ELEM_CLASSES  15
/* look-up table to translate the gmsh tree class to a t8code tree class.
 */
const t8_eclass_t   t8_msh_tree_type_to_eclass[T8_NUM_GMSH_ELEM_CLASSES + 1] = {
  T8_ECLASS_COUNT,              /* 0 is not valid */
  T8_ECLASS_LINE,               /* 1 */
  T8_ECLASS_TRIANGLE,
  T8_ECLASS_QUAD,
  T8_ECLASS_TET,
  T8_ECLASS_HEX,                /* 5 */
  T8_ECLASS_PRISM,
  T8_ECLASS_PYRAMID,            /* 7 This is the last first order tree type,
                                   except the Point, which is type 15 */
  /* We do not support type 8 to 14 */
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT,
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT,
  T8_ECLASS_VERTEX              /* 15 */
};

/* translate the msh file vertex number to the t8code vertex number.
 * See also http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering */
/* TODO: Check if these are correct */
const int           t8_msh_tree_vertex_to_t8_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 3, 2, 4, 5, 7, 6},     /* HEX */
  {0, 1, 2, 3},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* translate the t8code vertex number to the .msh file vertex number.
 * See also http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering */
/* TODO: Check if these are correct */
const int           t8_vertex_to_msh_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 3, 2, 4, 5, 7, 6},     /* HEX */
  {0, 1, 2, 3},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* TODO: if partitioned then only add the needed face-connections to join faces
 *       maybe also only trees and ghosts to classes.
 *       Specifying all face-connections makes commit algorithm slow! */

/* TODO: eventually compute neighbours only from .node and .ele files, since
 *       creating .neigh files with tetgen/triangle is not common and even seems
 *       to not work sometimes */

/* Read a the next line from a file stream that does not start with '#' or
 * contains only whitespaces (tabs etc.)
 *
 * \param [in,out] line     An allocated string to store the line.
 * \param [in,out] n        The number of allocated bytes.
 *                          If more bytes are needed line is reallocated and
 *                          the new number of bytes is stored in n.
 * \param [in]     fp       The file stream to read from.
 * \return                  The number of read arguments of the last line read.
 *                          negative on failure */
static int
t8_cmesh_msh_read_next_line (char **line, size_t *n, FILE *fp)
{
  int                 retval;

  do {
    /* read first non-comment line from file */
    /* TODO: getline depends on IEEE Std 1003.1-2008 (``POSIX.1'')
     *       p4est therefore has its own getline function in p4est_connectivity.h. */
    retval = getline (line, n, fp);
    if (retval < 0) {
      return retval;
    }
  }
  /* check if line is a comment (trailing '#') or consists solely of
   * blank spaces/tabs */
  while (*line[0] == '#' || strspn (*line, " \t\r\v\n") == strlen (*line));
  return retval;
}

/**
 * Reads an open msh-file and checks whether the MeshFormat-Version is supported by t8code or not. 
 * 
 * \param fp  A file pointer, pointing to an open gmsh-file
 * \return int if successfull, the version number is returned, an error-code if not.
 */
static int
t8_cmesh_check_version_of_msh_file (FILE *fp)
{
  char               *line = (char *) malloc (1024);
  char                first_word[2048] = "\0";
  size_t              linen = 1024;
  int                 retval;
  int                 version_number, sub_version_number;
  int                 check_format;
  int                 check_version = 0;

  T8_ASSERT (fp != NULL);

  /* Go to the beginning of the file. */
  fseek (fp, 0, SEEK_SET);

  /* Search for the line starting with "$MeshFormat". */
  while (!feof (fp) && strcmp (first_word, "$MeshFormat")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    retval = sscanf (line, "%2048s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf
        ("Reading the msh-file in order to check the MeshFormat-number failed.\n");
      goto die_format;
    }
  }

  /* Got to the next line containing the MeshFormat. */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Get the MeshFormat number of the file */
  retval =
    sscanf (line, "%d.%d %d", &version_number, &sub_version_number,
            &check_format);

  /*Checking for read/write error. */
  if (retval != 3) {
    t8_debugf ("Reading of the MeshFormat-number failed.\n");
    goto die_format;
  }

  /* Checks if the file is of Binary-type. */
  if (check_format) {
    t8_global_errorf
      ("Incompatible file-type. t8code works with ASCII-type msh-files with the versions:\n");
    for (int n_versions = 0;
         n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
      t8_global_errorf ("%d.X\n",
                        t8_cmesh_supported_msh_file_versions[n_versions]);
    }
    goto die_format;
  }

  /* Check if MeshFormat-number is compatible. */
  for (int n_versions = 0;
       n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
    if (version_number == t8_cmesh_supported_msh_file_versions[n_versions]) {
      check_version = 1;
    }
  }
  if (check_version) {
    t8_debugf ("This version of msh-file (%d.%d) is supported.\n",
               version_number, sub_version_number);
    free (line);
    return version_number;
  }
  else {
    t8_global_errorf
      ("This version of msh-file (%d.%d) is currently not supported by t8code, "
       "t8code supports ASCII files with the versions:\n",
       version_number, sub_version_number);
    for (int n_versions = 0;
         n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
      t8_global_errorf ("%d.X\n",
                        t8_cmesh_supported_msh_file_versions[n_versions]);
    }
    free (line);
    return 0;
  }

/* Will be executed, if reading the MeshFormat failed. */
die_format:
  /* Free memory. */
  free (line);
  /* Returning as error code. */
  return -1;
}

/* *INDENT-OFF* */
t8_read_status_t
t8_gmsh_reader::read (t8_cmesh_t cmesh)
/* *INDENT-ON* */

{
  t8_debugf ("[D] test read\n");
  return T8_READ_SUCCESS;
}

/* *INDENT-OFF* */
t8_read_status_t
t8_gmsh_reader::set_source (const t8_extern_t * source)
{
  if (source == NULL) {
    return T8_READ_FAIL;
  }
  else {
    char filepath[BUFSIZ];
    snprintf(filepath, BUFSIZ, "%s.msh", (const char *) source);
    /* Open the file */
    t8_debugf ("Opening file %s\n", filepath);
    file = fopen (filepath, "r");
    if(file == NULL){
        return T8_READ_FAIL;
    }
    /* Check if the msh-file version is compatible. */
    msh_version = t8_cmesh_check_version_of_msh_file (file);
    if(msh_version < 1 || (msh_version != 2 && msh_version != 4)){
        fclose(file);
        t8_debugf ("The reading process of the msh-file has failed and the file has been closed.\n");
        return T8_READ_FAIL;
    }
    return T8_READ_SUCCESS;
  }
}
/* *INDENT-ON* */

t8_gmsh_reader::t8_gmsh_reader ()
{
}

t8_gmsh_reader::~t8_gmsh_reader ()
{
}

#ifdef T8_ENABLE_DEBUG
int
t8_gmsh_reader::valid ()
{
  /* TODO: replace with something better as soon as more functionalitiy is implemented. */
  return 1;
}
#endif
