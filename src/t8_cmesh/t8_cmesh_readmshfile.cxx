/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8_eclass.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#if T8_ENABLE_OCC
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.h>
#endif
#include "t8_cmesh_types.h"
#include "t8_cmesh_stash.h"
#include <unordered_set>
#include <optional>
#include <vector>

#ifdef _WIN32
#include "t8_windows.h"
#endif

/* The supported number of gmesh tree classes.
 * Currently, we only support first order trees.
 */
#define T8_NUM_GMSH_ELEM_CLASSES 15
/* look-up table to translate the gmsh tree class to a t8code tree class.
 */
/* clang-format off */
const t8_eclass_t t8_msh_tree_type_to_eclass[T8_NUM_GMSH_ELEM_CLASSES + 1] = {
  T8_ECLASS_COUNT,     /* 0 is not valid */
  T8_ECLASS_LINE,      /* 1 */
  T8_ECLASS_TRIANGLE, 
  T8_ECLASS_QUAD, 
  T8_ECLASS_TET, 
  T8_ECLASS_HEX,        /* 5 */
  T8_ECLASS_PRISM, 
  T8_ECLASS_PYRAMID,    /* 7 This is the last first order tree type, except the Point, which is type 15 */
  /* We do not support type 8 to 14 */
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT, 
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT,
  T8_ECLASS_VERTEX      /* 15 */
};
/* clang-format on */

/* translate the msh file vertex number to the t8code vertex number.
 * See also http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering */
const int t8_msh_tree_vertex_to_t8_vertex_num[T8_ECLASS_COUNT][8] = {
  { 0 },                      /* VERTEX */
  { 0, 1 },                   /* LINE */
  { 0, 1, 3, 2 },             /* QUAD */
  { 0, 1, 2 },                /* TRIANGLE */
  { 0, 1, 3, 2, 4, 5, 7, 6 }, /* HEX */
  { 0, 1, 3, 2 },             /* TET */
  { 0, 1, 2, 3, 4, 5 },       /* PRISM */
  { 0, 1, 3, 2, 4 }           /* PYRAMID */
};

/* translate the t8code vertex number to the .msh file vertex number.
 * See also http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering */
const int t8_vertex_to_msh_vertex_num[T8_ECLASS_COUNT][8] = {
  { 0 },                      /* VERTEX */
  { 0, 1 },                   /* LINE */
  { 0, 1, 3, 2 },             /* QUAD */
  { 0, 1, 2 },                /* TRIANGLE */
  { 0, 1, 3, 2, 4, 5, 7, 6 }, /* HEX */
  { 0, 1, 3, 2 },             /* TET */
  { 0, 1, 2, 3, 4, 5 },       /* PRISM */
  { 0, 1, 3, 2, 4 }           /* PYRAMID */
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
  int retval;

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

/** The nodes are stored in the .msh file in the format
 *
 * $Nodes
 * n_nodes     // The number of nodes
 * i x_i y_i z_i  // the node index and the node coordinates
 * j x_j y_j z_j
 * .....
 * $EndNodes
 *
 * The node indices do not need to be in consecutive order.
 * We thus use a hash table to read all node indices and coordinates.
 * The hash value is the node index modulo the number of nodes.
 */
struct t8_msh_file_node
{
  /**
   * Default constructor.
   */
  t8_msh_file_node ()
    : parameters ({ -1, -1 }), coordinates ({ -1, -1 }), index (-1), parametric (false), entity_dim (-1),
      entity_tag (-1)
  {
  }

  /**
   * Constructor for non-parametric nodes.
   * \param [in, out] id        ID of the node.
   * \param [in, out] coords    Coords of the node.
   */
  t8_msh_file_node (t8_gloidx_t id, std::array<double, 3> coords)
    : parameters ({ -1, -1 }), coordinates (coords), index (id), parametric (false), entity_dim (-1), entity_tag (-1)
  {
  }

  /**
   * Constructor for parametric nodes.
   * \param [in, out] id        ID of the node.
   * \param [in, out] coords    Coords of the node.
   * \param [in, out] params    Parameters of the node in the parametric space.
   * \param [in] parametric True if the node is parametric, false otherwise.#
   * \param [in] entity_dim The dimension of the entity to which the node belongs.
   * \param [in] entity_tag The tag of the entity to which the node belongs.
   */
  t8_msh_file_node (t8_gloidx_t id, std::array<double, 3> coords, std::array<double, 2> params, bool parametric,
                    int entity_dim, t8_locidx_t entity_tag)
    : parameters (params), coordinates (coords), index (id), parametric (parametric), entity_dim (entity_dim),
      entity_tag (entity_tag)
  {
  }

  /**
   * Copy constructor
   * \param [in] other The node to copy.
   */
  t8_msh_file_node (const t8_msh_file_node &other)
    : parameters (other.parameters), coordinates (other.coordinates), index (other.index),
      parametric (other.parametric), entity_dim (other.entity_dim), entity_tag (other.entity_tag)
  {
  }

  /**
   * Move constructor.
   * \param [in] other The node to move.
   */
  t8_msh_file_node (t8_msh_file_node &&other) noexcept
    : parameters (std::move (other.parameters)), coordinates (std::move (other.coordinates)), index (other.index),
      parametric (other.parametric), entity_dim (other.entity_dim), entity_tag (other.entity_tag)
  {
    other.index = -1;
    other.parametric = false;
    other.entity_dim = -1;
    other.entity_tag = -1;
  }

  /**
   * Copy assignment operator.
   * \param [in] other The node to copy.
   * \return           Reference to this node.
   */
  t8_msh_file_node &
  operator= (const t8_msh_file_node &other)
  {
    if (this != &other) {
      parameters = other.parameters;
      coordinates = other.coordinates;
      index = other.index;
      parametric = other.parametric;
      entity_dim = other.entity_dim;
      entity_tag = other.entity_tag;
    }
    return *this;
  }

  /**
   * Move assignment operator.
   * \param [in] other The node to move.
   * \return           Reference to this node.
   */
  t8_msh_file_node &
  operator= (t8_msh_file_node &&other) noexcept
  {
    if (this != &other) {
      parameters = std::move (other.parameters);
      coordinates = std::move (other.coordinates);
      index = other.index;
      parametric = other.parametric;
      entity_dim = other.entity_dim;
      entity_tag = other.entity_tag;

      other.index = -1;
      other.parametric = false;
      other.entity_dim = -1;
      other.entity_tag = -1;
    }
    return *this;
  }

  std::array<double, 2> parameters;  /**< Parameters of the node in the parametric space, if applicable.
                                           * For example, for a point on a curve, this would be the parameter on the curve. */
  std::array<double, 3> coordinates; /**< Coordinates of the node in physical space. */
  t8_gloidx_t index;                 /**< The index of the node in the msh file. */
  bool parametric;                   /**< True if the node is parametric, false otherwise.
                                           * If true, the parameters are stored in the parameters array. */
  int entity_dim;                    /**< The dimension of the entity to which the node belongs.
                                           * For example, for a point on a curve, this would be 1. */
  t8_locidx_t entity_tag;            /**< The tag of the entity to which the node belongs.
                                           * For example, for a point on a curve, this would be the tag of the curve. */
};

/**
 * Hasher for msh file nodes.
 */
struct t8_msh_node_hasher
{
  /**
   * The number of nodes in the msh file.
   * This is used to compute the hash value.
   */
  t8_locidx_t num_nodes;

  /**
   * Constructor of the node hasher
   * \param [in, out] num_nodes   Number of nodes in the msh file.
   */
  t8_msh_node_hasher (t8_locidx_t num_nodes): num_nodes (num_nodes)
  {
  }

  /**
   * Hasher function.
   * \param [in, out] node        Msh node to create a hash value for.
   * \return                      The hash.
   */
  t8_gloidx_t
  operator() (const t8_msh_file_node &node) const
  {
    return node.index % num_nodes;
  }
};

/**
 * /struct t8_msh_node_equal
 *
 * Equality operator for msh file nodes.
 * This is used to compare nodes in the hash table.
 */
struct t8_msh_node_equal
{
  /**
   * Checks two nodes for equality.
   * \param [in, out] lhs       First node.
   * \param [in, out] rhs       Second node.
   * \return                    True if nodes are equal.
   */
  bool
  operator() (const t8_msh_file_node &lhs, const t8_msh_file_node &rhs) const
  {
    return lhs.index == rhs.index;
  }
};

/** Hashtable to store msh file nodes. */
typedef std::unordered_set<t8_msh_file_node, t8_msh_node_hasher, t8_msh_node_equal> t8_msh_node_table;

/** Vector which stores the vertex indices of each tree in the t8code order. */
typedef std::vector<std::vector<t8_gloidx_t>> t8_msh_tree_vertex_indices;

/* Reads an open msh-file and checks whether the MeshFormat-Version is supported by t8code or not. */
static int
t8_cmesh_check_version_of_msh_file (FILE *fp)
{
  char *line = (char *) malloc (1024);
  char first_word[2048] = "\0";
  size_t linen = 1024;
  int retval;
  int version_number, sub_version_number;
  int check_format;
  int check_version = 0;

  T8_ASSERT (fp != NULL);

  /* Go to the beginning of the file. */
  fseek (fp, 0, SEEK_SET);

  /* Search for the line starting with "$MeshFormat". */
  while (!feof (fp) && strcmp (first_word, "$MeshFormat")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    retval = sscanf (line, "%2047s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Reading the msh-file in order to check the MeshFormat-number failed.\n");
      goto die_format;
    }
  }

  /* Got to the next line containing the MeshFormat. */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Get the MeshFormat number of the file */
  retval = sscanf (line, "%d.%d %d", &version_number, &sub_version_number, &check_format);

  /*Checking for read/write error. */
  if (retval != 3) {
    t8_debugf ("Reading of the MeshFormat-number failed.\n");
    goto die_format;
  }

  /* Checks if the file is of Binary-type. */
  if (check_format) {
    t8_global_errorf ("Incompatible file-type. t8code works with ASCII-type msh-files with the versions:\n");
    for (int n_versions = 0; n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
      t8_global_errorf ("%d.X\n", t8_cmesh_supported_msh_file_versions[n_versions]);
    }
    goto die_format;
  }

  /* Check if MeshFormat-number is compatible. */
  for (int n_versions = 0; n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
    if (version_number == t8_cmesh_supported_msh_file_versions[n_versions]) {
      check_version = 1;
    }
  }
  if (check_version) {
    t8_debugf ("This version of msh-file (%d.%d) is supported.\n", version_number, sub_version_number);
    free (line);
    return version_number;
  }
  else {
    t8_global_errorf ("This version of msh-file (%d.%d) is currently not supported by t8code, "
                      "t8code supports ASCII files with the versions:\n",
                      version_number, sub_version_number);
    for (int n_versions = 0; n_versions < T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS; ++n_versions) {
      t8_global_errorf ("%d.X\n", t8_cmesh_supported_msh_file_versions[n_versions]);
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

/* Read an open .msh file of version 2 and parse the nodes into a hash table. */
static std::optional<t8_msh_node_table>
t8_msh_file_2_read_nodes (FILE *fp)
{
  t8_locidx_t ln, last_index;
  char *line = (char *) malloc (1024);
  char first_word[2048] = "\0";
  size_t linen = 1024;
  int retval;
  long lnum_nodes;

  T8_ASSERT (fp != NULL);
  /* Go to the beginning of the file */
  fseek (fp, 0, SEEK_SET);
  /* Search for the line beginning with "$Nodes" */
  while (!feof (fp) && strcmp (first_word, "$Nodes")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2047s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num nodes.\n");
      t8_debugf ("The line is %s.", line);
      free (line);
      return std::nullopt;
    }
  }
  /* Read the line containing the number of nodes */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Read the number of nodes in a long int before converting it
   * to t8_locidx_t. */
  retval = sscanf (line, "%li", &lnum_nodes);
  /* Checking for read/write error */
  if (retval != 1) {
    t8_global_errorf ("Premature end of line while reading num nodes.\n");
    t8_debugf ("The line is %s.", line);
    free (line);
    return std::nullopt;
  }

  /* Create the hash table */
  t8_msh_node_hasher hasher (lnum_nodes);
  t8_msh_node_table node_table (lnum_nodes, hasher);

  /* read each node and add it to the hash table */
  last_index = 0;
  std::array<double, 3> coords;
  t8_gloidx_t index;
  for (ln = 0; ln < lnum_nodes; ln++) {
    /* Read the next line. Its format should be %i %f %f %f
     * The node index followed by its coordinates. */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Error reading node file.\n");
      free (line);
      return std::nullopt;
    }
    /* Fill the node with the entries in the file */
    retval = sscanf (line, "%li %lf %lf %lf", &index, &coords[0], &coords[1], &coords[2]);
    if (retval != 4) {
      t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
      free (line);
      return std::nullopt;
    }
    /* Insert the node in the hash table */
    auto emplaced = node_table.emplace (t8_msh_file_node { index, coords });
    /* If second value is false then the node was already in the hash table.
     * This case should not occur. */
    if (emplaced.second == false) {
      t8_global_errorf ("Node %li defined more than once.\n", index);
      free (line);
      return std::nullopt;
    }
    T8_ASSERT (emplaced.first->index == index);
    last_index = index;
  }

  free (line);
  t8_debugf ("Successfully read all nodes.\n");
  return std::make_optional<t8_msh_node_table> (node_table);
}

/* Read an open .msh file of version 4 and parse the nodes into a hash table. */
static std::optional<t8_msh_node_table>
t8_msh_file_4_read_nodes (FILE *fp)
{
  t8_locidx_t ln, last_index;
  char *line = (char *) malloc (1024);
  char first_word[2048] = "\0";
  size_t linen = 1024;
  int retval;
  long num_nodes_in_block, lnum_nodes, lnum_blocks;
  long *index_buffer;

  T8_ASSERT (fp != NULL);
  /* Go to the beginning of the file */
  fseek (fp, 0, SEEK_SET);
  /* Search for the line beginning with "$Nodes" */
  while (!feof (fp) && strcmp (first_word, "$Nodes")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2047s", first_word);
    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading nodes.\n");
      t8_debugf ("The line is %s.", line);
      free (line);
      return std::nullopt;
    }
  }

  /* Read the line containing the number of entity block and nodes */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Read the number of nodes in a long int before converting it
   * to t8_locidx_t. */
  retval = sscanf (line, "%li %li %*i %*i", &lnum_blocks, &lnum_nodes);
  /* Checking for read/write error */
  if (retval != 2) {
    t8_global_errorf ("Premature end of line while reading num nodes and num blocks.\n");
    t8_debugf ("The line is %s.", line);
    free (line);
    return std::nullopt;
  }

  /* Create the hash table */
  t8_msh_node_hasher hasher (lnum_nodes);
  t8_msh_node_table node_table (lnum_nodes, hasher);

  /* read each node and add it to the hash table */
  last_index = 0;
  std::array<double, 3> coords;
  std::array<double, 2> params;
  int entity_dim;
  t8_locidx_t entity_tag;
  int parametric;
  for (long n_block = 0; n_block < lnum_blocks; ++n_block) {
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Error reading node file.\n");
      free (line);
      return std::nullopt;
    }
    retval = sscanf (line, "%i %i %i %li", &entity_dim, &entity_tag, &parametric, &num_nodes_in_block);
    if (retval != 4) {
      t8_global_errorf ("Error reading block after node %li in $Nodes section.\n", (long) last_index);
      free (line);
      return std::nullopt;
    }
    /* Read all node indices in this block */
    index_buffer = T8_ALLOC (long, num_nodes_in_block);
    for (ln = 0; ln < num_nodes_in_block; ++ln) {
      retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
      if (retval < 0) {
        t8_global_errorf ("Error reading node file.\n");
        free (line);
        return std::nullopt;
      }
      retval = sscanf (line, "%li", &index_buffer[ln]);
      if (retval != 1) {
        t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
        free (line);
        return std::nullopt;
      }
    }
    /* Read all coordinates and parameters in this block */
    for (ln = 0; ln < num_nodes_in_block; ++ln) {
      /* Read the next line. Its format should be
       * %f %f %f         if not parameterized,
       * %f %f %f %f      if parameterized and entity_dim == 1,
       * %f %f %f %f %f   if parameterized and entity_dim == 2.
       * The coordinates followed by their parameters. */
      retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
      if (retval < 0) {
        t8_global_errorf ("Error reading node file\n");
        free (line);
        return std::nullopt;
      }
      /* Fill the node with the entries in the file */
      if (!parametric) {
        retval = sscanf (line, "%lf %lf %lf", &coords[0], &coords[1], &coords[2]);
        if (retval != 3) {
          t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
          free (line);
          return std::nullopt;
        }
      }
      else {
        /* Check for entity_dim and retrieve parameters accordingly */
        switch (entity_dim) {
        case 1:
          retval = sscanf (line, "%lf %lf %lf %lf", &coords[0], &coords[1], &coords[2], &params[0]);
          if (retval != 4) {
            t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
            free (line);
            return std::nullopt;
          }
          break;
        case 2:
          retval = sscanf (line, "%lf %lf %lf %lf %lf", &coords[0], &coords[1], &coords[2], &params[0], &params[1]);
          if (retval != 5) {
            t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
            free (line);
            return std::nullopt;
          }
          break;
        default:
          t8_global_errorf ("Error reading node file after node %li.\n", (long) last_index);
          free (line);
          return std::nullopt;
        }
      }
      /* Insert the node in the hash table */
      auto emplaced = node_table.emplace (
        t8_msh_file_node { index_buffer[ln], coords, params, (bool) parametric, entity_dim, entity_tag });
      /* If second value is false then the node was already in the hash table.
     * This case should not occur. */
      if (emplaced.second == false) {
        t8_global_errorf ("Node %li defined more than once.\n", index_buffer[ln]);
        free (line);
        return std::nullopt;
      }
      T8_ASSERT (emplaced.first->index == index_buffer[ln]);
      last_index = index_buffer[ln];
    }
    T8_FREE (index_buffer);
  }
  free (line);
  t8_debugf ("Successfully read all Nodes.\n");
  return std::make_optional<t8_msh_node_table> (node_table);
}

/**
 * Adds the elements of \a fp and dimension \a dim into the \a cmesh.
 * Returns a list of all vertex indices of each tree. 
 * \param [in, out] cmesh     The cmesh.
 * \param [in, out] fp        The msh file.
 * \param [in, out] vertices  A hashtable filled with the nodes of the msh file.
 * \param [in, out] dim       The dimension of nodes to read in.
 * \return 
 */
static std::optional<t8_msh_tree_vertex_indices>
t8_cmesh_msh_file_2_read_eles (t8_cmesh_t cmesh, FILE *fp, const t8_msh_node_table vertices, const int dim)
{
  char *line = (char *) malloc (1024), *line_modify;
  char first_word[2048] = "\0";
  size_t linen = 1024;
  t8_locidx_t num_trees, tree_loop;
  t8_gloidx_t tree_count;
  t8_eclass_t eclass;
  t8_msh_file_node Node;
  long lnum_trees;
  int retval;
  int ele_type, num_tags;
  int num_nodes;
  std::array<double, T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM> tree_vertices;

  T8_ASSERT (fp != NULL);
  /* Search for the line beginning with "$Elements" */
  while (!feof (fp) && strcmp (first_word, "$Elements")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2047s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num trees.\n");
      t8_debugf ("The line is %s", line);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
  }

  /* Read the line containing the number of trees */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Since t8_locidx_t could be int32 or int64, we first read the
   * number of trees in a long int and store it as t8_locidx_t later. */
  retval = sscanf (line, "%li", &lnum_trees);
  /* Checking for read/write error */
  if (retval != 1) {
    t8_global_errorf ("Premature end of line while reading num trees.\n");
    t8_debugf ("The line is %s", line);
    free (line);
    t8_cmesh_destroy (&cmesh);
    return std::nullopt;
  }
  num_trees = lnum_trees;
  /* Check for type conversion error */
  T8_ASSERT (num_trees == lnum_trees);

  /* Reserve memory for vertex indices */
  t8_msh_tree_vertex_indices vertex_indices (num_trees);

  tree_count = 0; /* The index of the next tree to insert */
  for (tree_loop = 0; tree_loop < num_trees; tree_loop++) {
    /* Read the next line containing tree information */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Premature end of line while reading trees.\n");
      t8_debugf ("The line is %s", line);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
    /* The line describing the tree looks like
     * tree_number tree_type Number_tags tag_1 ... tag_n Node_1 ... Node_m
     *
     * We ignore the tree number, read the type and the number of (integer) tags.
     * We also ignore the tags and after we know the type, we read the
     * nodes.
     */
    sscanf (line, "%*i %i %i", &ele_type, &num_tags);
    /* Check if the tree type is supported */
    if (ele_type > T8_NUM_GMSH_ELEM_CLASSES || ele_type < 0
        || t8_msh_tree_type_to_eclass[ele_type] == T8_ECLASS_COUNT) {
      t8_global_errorf ("Tree type %i is not supported by t8code.\n", ele_type);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
    /* Continue if tree type is supported */
    eclass = t8_msh_tree_type_to_eclass[ele_type];
    T8_ASSERT (eclass != T8_ECLASS_COUNT);

    if (t8_eclass_to_dimension[eclass] > dim) {
      t8_debugf (
        "Warning: Encountered an element with a dimension higher than %d. Did you set the correct dimension?\n", dim);
    }

    /* Check if the tree is of the correct dimension */
    if (t8_eclass_to_dimension[eclass] == dim) {
      /* The tree is of the correct dimension,
       * add it to the cmesh and read its nodes */
      t8_cmesh_set_tree_class (cmesh, tree_count, eclass);
      line_modify = line;
      /* Since the tags are stored before the node indices, we need to
       * skip them first. But since the number of them is unknown and the
       * length (in characters) of them, we have to skip one by one. */
      for (int i_tag = 0; i_tag < 3 + num_tags; i_tag++) {
        T8_ASSERT (strcmp (line_modify, "\0"));
        /* move line_modify to the next word in the line */
        (void) strsep (&line_modify, " ");
      }
      /* At this point line_modify contains only the node indices. */
      num_nodes = t8_eclass_num_vertices[eclass];
      std::vector<t8_gloidx_t> node_indices (num_nodes, -1);
      for (int i_node = 0; i_node < num_nodes; i_node++) {
        const int t8_vertex_num = t8_msh_tree_vertex_to_t8_vertex_num[eclass][i_node];
        T8_ASSERT (strcmp (line_modify, "\0"));
        retval = sscanf (line_modify, "%li", &node_indices[t8_vertex_num]);
        if (retval != 1) {
          t8_global_errorf ("Premature end of line while reading tree.\n");
          t8_debugf ("The line is %s", line);
          free (line);
          t8_cmesh_destroy (&cmesh);
          return std::nullopt;
        }

        /* Get node from the hashtable */
        Node.index = node_indices[t8_vertex_num];
        const auto found_node = vertices.find (Node);
        if (found_node == vertices.end ()) {
          t8_global_errorf ("Could not find Node %li.\n", node_indices[t8_vertex_num]);
          free (line);
          t8_cmesh_destroy (&cmesh);
          return std::nullopt;
        }

        /* Add node coordinates to the tree vertices */
        tree_vertices[3 * t8_vertex_num] = found_node->coordinates[0];
        tree_vertices[3 * t8_vertex_num + 1] = found_node->coordinates[1];
        tree_vertices[3 * t8_vertex_num + 2] = found_node->coordinates[2];
        /* move line_modify to the next word in the line */
        (void) strsep (&line_modify, " ");
      }

      /* Add the node indices to return vector. */
      vertex_indices[tree_count] = std::move (node_indices);
      /* Detect and correct negative volumes */
      if (t8_cmesh_tree_vertices_negative_volume (eclass, tree_vertices.data (), num_nodes)) {
        /* The volume described is negative. We need to change vertices.
         * For tets we switch 0 and 3.
         * For prisms we switch 0 and 3, 1 and 4, 2 and 5.
         * For hexahedra we switch 0 and 4, 1 and 5, 2 and 6, 3 and 7.
         * For pyramids we switch 0 and 4 */
        double temp;
        int num_switches = 0;
        int switch_indices[4] = { 0 };
        int iswitch;
        T8_ASSERT (t8_eclass_to_dimension[eclass] > 1);
        t8_debugf ("Correcting negative volume of tree %li\n", static_cast<long> (tree_count));
        switch (eclass) {
        case T8_ECLASS_TRIANGLE:
        case T8_ECLASS_QUAD:
          /* We switch vertex 1 and vertex 2. */
          num_switches = 2;
          switch_indices[0] = 0;
          switch_indices[1] = 2;
          break;
        case T8_ECLASS_TET:
          /* We switch vertex 0 and vertex 3 */
          num_switches = 1;
          switch_indices[0] = 3;
          break;
        case T8_ECLASS_PRISM:
          num_switches = 3;
          switch_indices[0] = 3;
          switch_indices[1] = 4;
          switch_indices[2] = 5;
          break;
        case T8_ECLASS_HEX:
          num_switches = 4;
          switch_indices[0] = 4;
          switch_indices[1] = 5;
          switch_indices[2] = 6;
          switch_indices[3] = 7;
          break;
        case T8_ECLASS_PYRAMID:
          num_switches = 1;
          switch_indices[0] = 4;
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }

        for (iswitch = 0; iswitch < num_switches; ++iswitch) {
          /* We switch vertex 0 + iswitch and vertex switch_indices[iswitch] */
          for (int i_dim = 0; i_dim < T8_ECLASS_MAX_DIM; i_dim++) {
            temp = tree_vertices[3 * iswitch + i_dim];
            tree_vertices[3 * iswitch + i_dim] = tree_vertices[3 * switch_indices[iswitch] + i_dim];
            tree_vertices[3 * switch_indices[iswitch] + i_dim] = temp;
          }
        }
        T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume (eclass, tree_vertices.data (), num_nodes));
      } /* End of negative volume handling */
      /* Set the vertices of this tree */
      t8_cmesh_set_tree_vertices (cmesh, tree_count, tree_vertices.data (), num_nodes);
    }
    /* advance the tree counter */
    tree_count++;
  }
  free (line);
  if (tree_count == 0) {
    t8_global_errorf ("Warning: No %iD elements found in msh file.\n", dim);
    t8_cmesh_destroy (&cmesh);
    return std::nullopt;
  }
  return std::make_optional<t8_msh_tree_vertex_indices> (vertex_indices);
}

#if T8_ENABLE_OCC
/** Corrects the parameters on closed geometries to prevent disorted elements.
 * \param [in]      geometry_dim    The dimension of the geometry.
 *                                  1 for edges, 2 for surfaces.
 * \param [in]      geometry_index  The index of the geometry.
 * \param [in]      num_face_nodes  The number of the nodes of the surface.
 *                                  NULL if the geometry is an edge.
 * \param [in]      geometry_cad    The cad_geometry.
 * \param [in,out]  parameters      The parameters to be corrected.
 */
static void
t8_cmesh_correct_parameters_on_closed_geometry (const int geometry_dim, const int geometry_index,
                                                const int num_face_nodes, const t8_geometry_cad *geometry_cad,
                                                double *parameters)
{
  switch (geometry_dim) {
    /* Check for closed U parameter in case of an edge. */
  case 1:
    /* Only correct the U parameter if the edge is closed. */
    if (geometry_cad->get_cad_manager ()->t8_geom_is_edge_closed (geometry_index)) {
      /* Get the parametric bounds of the closed geometry 
       * edge    -> [Umin, Umax]
       */
      double parametric_bounds[2];
      /* Get the parametric edge bounds. */
      geometry_cad->get_cad_manager ()->t8_geom_get_edge_parametric_bounds (geometry_index, parametric_bounds);
      /* Check the upper an the lower parametric bound. */
      for (int bound = 0; bound < 2; ++bound) {
        /* Iterate over both nodes of the edge. */
        for (int i_nodes = 0; i_nodes < 2; ++i_nodes) {
          /* Check if one of the U parameters lies on one of the parametric bounds. */
          if (std::abs (parameters[i_nodes] - parametric_bounds[bound]) <= T8_PRECISION_EPS) {
            /* Check the U parameter of the other node ((i_node + 1) % 2) to find out
             * to which parametric bound the tree is closer.
             */
            if (std::abs (parameters[(i_nodes + 1) % 2] - parametric_bounds[bound]) > T8_PRECISION_EPS) {
              /* Now check if the difference of the parameters of both nodes are bigger than the half parametric range.
               * In this case, the parameter at i_nodes has to be changed to the other parametric bound ((bound + 1) % 2).
               */
              if (std::abs (parameters[(i_nodes + 1) % 2] - parameters[i_nodes])
                  > ((parametric_bounds[1] - parametric_bounds[0]) / 2)) {
                /* Switch to the other parametric bound. */
                parameters[i_nodes] = parametric_bounds[(bound + 1) % 2];
                break;
              }
            }
          }
        }
      }
    }
    break;

    /* Check for closed U parameter and closed V parameter in case of a surface. */
  case 2:
    /* Iterate over both parameters. 0 stands for the U parameter an 1 for the V parameter. */
    for (int param_dim = 0; param_dim < 2; ++param_dim) {
      /* Only correct the surface parameters if they are closed */
      if (geometry_cad->get_cad_manager ()->t8_geom_is_surface_closed (geometry_index, param_dim)) {
        /* Get the parametric bounds of the closed geometry
         * surface -> [Umin, Umax, Vmin, Vmax]
         */
        double parametric_bounds[4];
        geometry_cad->get_cad_manager ()->t8_geom_get_face_parametric_bounds (geometry_index, parametric_bounds);
        /* Check the upper an the lower parametric bound. */
        for (int bound = 0; bound < 2; ++bound) {
          /* Iterate over every corner node of the tree. */
          for (int i_nodes = 0; i_nodes < num_face_nodes; ++i_nodes) {
            /* Check if one of the U parameters lies on one of the parametric bounds. */
            if (std::abs (parameters[2 * i_nodes + param_dim] - parametric_bounds[bound + 2 * param_dim])
                <= T8_PRECISION_EPS) {
              /* Iterate over every corner node of the tree again. */
              for (int j_nodes = 0; j_nodes < num_face_nodes; ++j_nodes) {
                /* Search for a U parameter that is non of the parametric bounds. To check
                 * whether the tree is closer to the lower or the upper parametric bound.
                 */
                if (std::abs (parameters[2 * j_nodes + param_dim] - parametric_bounds[bound + 2 * param_dim])
                      > T8_PRECISION_EPS
                    && std::abs (parameters[2 * j_nodes + param_dim]
                                 - parametric_bounds[((bound + 1) % 2) + 2 * param_dim])
                         > T8_PRECISION_EPS) {
                  /* Now check if the difference of the parameters of both nodes are bigger than the half parametric range.
                   * In this case, the parameter at i_nodes has to be changed to the other parametric bound ((bound + 1) % 2).
                   */
                  if (std::abs (parameters[2 * j_nodes + param_dim] - parameters[2 * i_nodes + param_dim])
                      > ((parametric_bounds[1 + 2 * param_dim] - parametric_bounds[0 + 2 * param_dim]) / 2)) {
                    /* Switch to the other parametric bound. */
                    parameters[2 * i_nodes + param_dim] = parametric_bounds[((bound + 1) % 2) + 2 * param_dim];
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
}
#endif /* T8_ENABLE_OCC */

/**
 * This function stores the entity dimensions, tags, and parameters of each node in a given 
 * tree in arrays. These arrays are then passed to the tree as attributes in cmesh.
 * \param [in, out] cmesh        The cmesh.
 * \param [in] tree_count   The index of the tree.
 * \param [in] tree_nodes   The array of the nodes of the tree.
 * \param [in] num_nodes    The number of nodes. 
 */
static void
t8_store_element_node_data (t8_cmesh_t cmesh, t8_gloidx_t tree_count,
                            std::array<t8_msh_file_node, T8_ECLASS_MAX_CORNERS> *tree_nodes, int num_nodes)
{

  int entity_dim_array[T8_ECLASS_MAX_CORNERS * 2];
  double parameter_array[T8_ECLASS_MAX_CORNERS * 2];

  for (int node_i = 0; node_i < num_nodes; node_i++) {
    t8_msh_file_node &current_node = (*tree_nodes)[node_i];

    /* Store entity_dim and the entity_tag in the entity_dim_array. */
    entity_dim_array[2 * node_i] = current_node.entity_dim;
    entity_dim_array[2 * node_i + 1] = current_node.entity_tag;

    /* Store the parameters of the node in the parameter_array. */
    parameter_array[2 * node_i] = current_node.parameters[0];
    parameter_array[2 * node_i + 1] = current_node.parameters[1];
  }

  /* Give the tree information about the dimension, tag and the parameters of the vertices. 
          * Each parameter set is given to the tree via its attribute key*/
  t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (), T8_CMESH_NODE_GEOMETRY_ATTRIBUTE_KEY,
                          entity_dim_array, T8_ECLASS_MAX_CORNERS * 2 * sizeof (int), 0);

  t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (), T8_CMESH_NODE_PARAMETERS_ATTRIBUTE_KEY,
                          parameter_array, T8_ECLASS_MAX_CORNERS * 2 * sizeof (int), 0);
}

#if T8_ENABLE_OCC
/**
 * This function calculates the parametric geometries of a tree based on its element class
 * and links it to either a CAD-based geometry or a linear geometry. It validates the element class
 * and assigns geometric attributes (e.g., face and edge parameters) to the tree. If the geometry
 * cannot be processed, the function returns 0.
 *
 * \param [in, out] cmesh The computational mesh to which the tree belongs.
 * \param [in] eclass The element class of the tree (e.g., triangle, quadrilateral, tetrahedron).
 * \param [in] dim The dimension of the geometry (e.g., 2D or 3D).
 * \param [in] tree_count The index of the tree in the computational mesh.
 * \param [in] cad_geometry_base A pointer to the CAD-based geometry object.
 * \param [in] linear_geometry_base A pointer to the linear geometry object.
 * \param [in] tree_nodes An array of nodes representing the vertices of the tree.
 * \param [in] face_nodes An array of nodes representing the faces of the tree.
 * \param [in] edge_nodes An array of nodes representing the edges of the tree.
 *
 * \return std::optional<t8_msh_tree_vertex_indices>
 *         - Returns a valid `t8_msh_tree_vertex_indices` object if the tree geometry is successfully processed.
 *         - Returns 0 if the geometry processing fails.
 */
static bool
t8_cmesh_process_tree_geometry (t8_cmesh_t cmesh, t8_eclass_t eclass, int dim, t8_gloidx_t tree_count,
                                const t8_geometry_c *cad_geometry_base, const t8_geometry_c *linear_geometry_base,
                                std::array<t8_msh_file_node, T8_ECLASS_MAX_CORNERS> tree_nodes,
                                std::array<t8_msh_file_node, T8_ECLASS_MAX_CORNERS_2D> face_nodes,
                                std::array<t8_msh_file_node, 2> edge_nodes)
{
  /* Calculate the parametric geometries of the tree */
  T8_ASSERT (cad_geometry_base->t8_geom_get_type () == T8_GEOMETRY_TYPE_CAD);
  const t8_geometry_cad *cad_geometry = dynamic_cast<const t8_geometry_cad *> (cad_geometry_base);
  /* Check for right element class */
  if (eclass != T8_ECLASS_TRIANGLE && eclass != T8_ECLASS_QUAD && eclass != T8_ECLASS_HEX && eclass != T8_ECLASS_TET
      && eclass != T8_ECLASS_PRISM) {
    t8_errorf ("%s element detected. The cad geometry currently only supports quad, tri, hex, tet and prism elements.",
               t8_eclass_to_string[eclass]);
    return 0;
  }
  int tree_is_linked = 0;
  double parameters[T8_ECLASS_MAX_CORNERS_2D * 2];
  int edge_geometries[T8_ECLASS_MAX_EDGES * 2] = { 0 };
  int face_geometries[T8_ECLASS_MAX_FACES] = { 0 };
  /* We look at each face to check, if it is linked to a cad surface */
  T8_ASSERT (t8_eclass_to_dimension[eclass] == dim);
  int num_faces;
  switch (dim) {
  case 0:
    num_faces = 0;
    break;
  case 1:
    num_faces = 0;
    break;
  case 2:
    num_faces = 1;
    break;
  case 3:
    num_faces = t8_eclass_num_faces[eclass];
    break;
  default:
    SC_ABORTF ("Invalid dimension of tree. Dimension: %i\n", dim);
  }
  for (int i_tree_faces = 0; i_tree_faces < num_faces; ++i_tree_faces) {
    const int face_eclass = dim == 2 ? eclass : t8_eclass_face_types[eclass][i_tree_faces];
    const int num_face_nodes = t8_eclass_num_vertices[face_eclass];
    const int num_face_edges = t8_eclass_num_faces[face_eclass];

    /* Save each node of face separately. Face nodes of 2D elements are also tree nodes.
             * Face nodes of 3D elements need to be translated to tree nodes. */
    for (int i_face_node = 0; i_face_node < num_face_nodes; ++i_face_node) {
      if (dim == 2) {
        face_nodes[i_face_node] = tree_nodes[i_face_node];
      }
      else {
        face_nodes[i_face_node] = tree_nodes[t8_face_vertex_to_tree_vertex[eclass][i_tree_faces][i_face_node]];
      }
    }

    /* A face can only be linked to an cad surface if all nodes of the face are parametric or on a vertex 
             * (gmsh labels nodes on vertices as not parametric) */
    int all_parametric = 1;
    for (int i_face_nodes = 0; i_face_nodes < num_face_nodes; ++i_face_nodes) {
      if (!face_nodes[i_face_nodes].parametric && face_nodes[i_face_nodes].entity_dim != 0) {
        all_parametric = 0;
        break;
      }
    }
    /* Skip face if not all nodes are parametric */
    if (!all_parametric) {
      continue;
    }
    /* Now we can check if the face is connected to a surface */
    int surface_index = 0;
    /* If one node is already on a surface we can check if the rest lies also on the surface. */
    for (int i_face_nodes = 0; i_face_nodes < num_face_nodes; ++i_face_nodes) {
      if (face_nodes[i_face_nodes].entity_dim == 2) {
        surface_index = face_nodes[i_face_nodes].entity_tag;
        break;
      }
    }
    /* If not we can take two curves and look if they share a surface and then use this surface */
    if (!surface_index) {
      /* To do this we can look if there are two curves, otherwise we have to check which vertices 
               * share the same curve. */
      int edge1_index = 0;
      int edge2_index = 0;
      /* We search for 2 different curves */
      for (int i_face_nodes = 0; i_face_nodes < num_face_nodes; ++i_face_nodes) {
        if (face_nodes[i_face_nodes].entity_dim == 1) {
          if (edge1_index == 0) {
            edge1_index = face_nodes[i_face_nodes].entity_tag;
          }
          else if (face_nodes[i_face_nodes].entity_tag != edge1_index) {
            edge2_index = face_nodes[i_face_nodes].entity_tag;
            break;
          }
        }
      }
      /* If there are less than 2 curves we can look at the vertices and check, 
               * if two of them are on the same curve */
      if (edge2_index == 0) {
        /* For each edge of face */
        for (int i_face_edges = 0; i_face_edges < num_face_edges; ++i_face_edges) {
          /* Save nodes separately */
          const int node1_number = t8_face_vertex_to_tree_vertex[face_eclass][i_tree_faces][0];
          const t8_msh_file_node node1 = face_nodes[node1_number];
          const int node2_number = t8_face_vertex_to_tree_vertex[face_eclass][i_tree_faces][1];
          const t8_msh_file_node node2 = face_nodes[node2_number];

          /* If both nodes are on a vertex we look if both vertices share an edge */
          if (node1.entity_dim == 0 && node2.entity_dim == 0) {
            int common_edge
              = cad_geometry->get_cad_manager ()->t8_geom_get_common_edge (node1.entity_tag, node2.entity_tag);
            if (common_edge > 0) {
              if (edge1_index == 0) {
                edge1_index = common_edge;
              }
              else if (edge2_index == 0 && common_edge != edge1_index) {
                edge2_index = common_edge;
                break;
              }
            }
          }
        }
      }
      if (edge2_index > 0) {
        surface_index = cad_geometry->get_cad_manager ()->t8_geom_get_common_face (edge1_index, edge2_index);
      }
      else {
        continue;
      }
    }
    /* Now we can check if every node lies on the surface and retrieve its parameters */
    if (surface_index) {
      int all_nodes_on_surface = 1;
      for (int i_face_nodes = 0; i_face_nodes < num_face_nodes; ++i_face_nodes) {
        /* We check if the node is on the right surface */
        if (face_nodes[i_face_nodes].entity_dim == 2) {
          /* Check if node is on the right surface */
          if (face_nodes[i_face_nodes].entity_tag != surface_index) {
            all_nodes_on_surface = 0;
            break;
          }
        }
        else {
          /* If it is on another geometry we retrieve its parameters */
          if (face_nodes[i_face_nodes].entity_dim == 0) {
            if (cad_geometry->get_cad_manager ()->t8_geom_is_vertex_on_face (face_nodes[i_face_nodes].entity_tag,
                                                                             surface_index)) {
              cad_geometry->get_cad_manager ()->t8_geom_get_parameters_of_vertex_on_face (
                face_nodes[i_face_nodes].entity_tag, surface_index, face_nodes[i_face_nodes].parameters.data ());
              face_nodes[i_face_nodes].entity_dim = 2;
            }
            else {
              all_nodes_on_surface = 0;
              break;
            }
          }
          if (face_nodes[i_face_nodes].entity_dim == 1) {
            if (cad_geometry->get_cad_manager ()->t8_geom_is_edge_on_face (face_nodes[i_face_nodes].entity_tag,
                                                                           surface_index)) {
              cad_geometry->get_cad_manager ()->t8_geom_edge_parameter_to_face_parameters (
                face_nodes[i_face_nodes].entity_tag, surface_index, num_face_nodes,
                face_nodes[i_face_nodes].parameters[0], NULL, face_nodes[i_face_nodes].parameters.data ());
              face_nodes[i_face_nodes].entity_dim = 2;
            }
            else {
              all_nodes_on_surface = 0;
              break;
            }
          }
        }
      }
      /* Abort if not all nodes are on the surface or if the surface is a plane */
      if (!all_nodes_on_surface || cad_geometry->get_cad_manager ()->t8_geom_is_plane (surface_index)) {
        continue;
      }
      /* If we have found a surface we link it to the face */
      face_geometries[i_tree_faces] = surface_index;
      tree_is_linked = 1;
      for (int i_face_edges = 0; i_face_edges < num_face_edges; ++i_face_edges) {
        /* We lock the edges of the face for surfaces, so that we do not link the same surface again 
                 * to the edges of the face */
        if (dim == 2) /* 2D */
        {
          edge_geometries[i_face_edges + t8_eclass_num_edges[eclass]] = -1;
        }
        else /* 3D */
        {
          edge_geometries[t8_face_edge_to_tree_edge[eclass][i_tree_faces][i_face_edges] + t8_eclass_num_edges[eclass]]
            = -1;
        }
      }
      /* We retrieve the parameters of the nodes and give them to the tree */
      for (int i_face_nodes = 0; i_face_nodes < num_face_nodes; ++i_face_nodes) {
        parameters[i_face_nodes * 2] = face_nodes[i_face_nodes].parameters[0];
        parameters[i_face_nodes * 2 + 1] = face_nodes[i_face_nodes].parameters[1];
      }
      /* Corrects the parameters on the surface if it is closed to prevent disorted elements. */
      for (int param_dim = 0; param_dim < 2; ++param_dim) {
        if (cad_geometry->get_cad_manager ()->t8_geom_is_surface_closed (surface_index, param_dim)) {
          t8_cmesh_correct_parameters_on_closed_geometry (2, surface_index, num_face_nodes, cad_geometry, parameters);
        }
      }

      t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (),
                              T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + i_tree_faces, parameters,
                              num_face_nodes * 2 * sizeof (double), 0);
    }
  }
  const int num_edges = t8_eclass_num_edges[eclass];
  /* Then we look for geometries linked to the edges */
  for (int i_tree_edges = 0; i_tree_edges < num_edges; ++i_tree_edges) {
    if (t8_eclass_to_dimension[eclass] == 3) {
      edge_nodes[0] = tree_nodes[t8_edge_vertex_to_tree_vertex[eclass][i_tree_edges][0]];
      edge_nodes[1] = tree_nodes[t8_edge_vertex_to_tree_vertex[eclass][i_tree_edges][1]];
    }
    else {
      edge_nodes[0] = tree_nodes[t8_face_vertex_to_tree_vertex[eclass][i_tree_edges][0]];
      edge_nodes[1] = tree_nodes[t8_face_vertex_to_tree_vertex[eclass][i_tree_edges][1]];
    }
    /* Both nodes have to be parametric or on a vertex to be linked to a curve or surface */
    if ((!edge_nodes[0].parametric && edge_nodes[0].entity_dim != 0)
        || (!edge_nodes[1].parametric && edge_nodes[1].entity_dim != 0)) {
      continue;
    }
    /* An edge can be linked to a curve as well as a surface. 
             * Therefore, we have to save the geometry dim and tag */
    int edge_geometry_dim = 0;
    int edge_geometry_tag = 0;
    /* We check which is the highest dim a node geometry has and what is its tag */
    if (edge_nodes[0].entity_dim > edge_nodes[1].entity_dim) {
      edge_geometry_dim = edge_nodes[0].entity_dim;
      if (edge_nodes[0].entity_dim > 0) {
        edge_geometry_tag = edge_nodes[0].entity_tag;
      }
    }
    else {
      edge_geometry_dim = edge_nodes[1].entity_dim;
      if (edge_nodes[1].entity_dim > 0) {
        edge_geometry_tag = edge_nodes[1].entity_tag;
      }
    }
    /* If both nodes are on two different faces we can skip this edge. */
    if (edge_nodes[0].entity_dim == 2 && edge_nodes[1].entity_dim == 2
        && edge_nodes[0].entity_tag != edge_nodes[1].entity_tag) {
      continue;
    }

    /* If one vertex lies on a geometry of a higher dim as the other, we have to check,
             * if the geometry of lower dimension is on that geometry. */
    {
      int is_on_geom = 1;
      for (int i_edge = 0; i_edge < 2; ++i_edge) {
        if (edge_geometry_dim == 2 && edge_nodes[i_edge].entity_dim == 1) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_edge_on_face (edge_nodes[i_edge].entity_tag,
                                                                          edge_geometry_tag)) {
            is_on_geom = 0;
            break;
          }
        }
        else if (edge_geometry_dim == 2 && edge_nodes[i_edge].entity_dim == 0) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_vertex_on_face (edge_nodes[i_edge].entity_tag,
                                                                            edge_geometry_tag)) {
            is_on_geom = 0;
            break;
          }
        }
        else if (edge_geometry_dim == 1 && edge_nodes[i_edge].entity_dim == 0) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_vertex_on_edge (edge_nodes[i_edge].entity_tag,
                                                                            edge_geometry_tag)) {
            is_on_geom = 0;
            break;
          }
        }
      }
      if (!is_on_geom) {
        continue;
      }
    }

    /* If both nodes are on a vertex we still got no edge. 
             * But we can look if both vertices share an edge and use this edge. 
             * If not we can skip this edge. */
    if (edge_geometry_dim == 0 && edge_geometry_tag == 0) {
      int common_curve = cad_geometry->get_cad_manager ()->t8_geom_get_common_edge (edge_nodes[0].entity_tag,
                                                                                    edge_nodes[1].entity_tag);
      if (common_curve > 0) {
        edge_geometry_tag = common_curve;
        edge_geometry_dim = 1;
      }
      else {
        continue;
      }
    }
    /* If both nodes are on different edges we have to look if both edges share a surface. 
             * If not we can skip this edge */
    if (edge_nodes[0].entity_dim == 1 && edge_nodes[1].entity_dim == 1
        && edge_nodes[0].entity_tag != edge_nodes[1].entity_tag) {
      int common_surface = cad_geometry->get_cad_manager ()->t8_geom_get_common_face (edge_nodes[0].entity_tag,
                                                                                      edge_nodes[1].entity_tag);
      if (common_surface > 0) {
        edge_geometry_tag = common_surface;
        edge_geometry_dim = 2;
      }
      else {
        continue;
      }
    }
    /* If we have found a curve we can look for the parameters */
    if (edge_geometry_dim == 1) {
      /* Check if adjacent faces carry a surface and if this edge lies on the surface */
      for (int i_adjacent_face = 0; i_adjacent_face < 2; ++i_adjacent_face) {
        if (face_geometries[t8_edge_to_face[eclass][i_tree_edges][i_adjacent_face]] > 0) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_edge_on_face (
                edge_geometry_tag, face_geometries[t8_edge_to_face[eclass][i_tree_edges][i_adjacent_face]])) {
            t8_global_errorf ("Error: Adjacent edge and face of a tree carry "
                              "incompatible geometries.\n");
            return 0;
          }
        }
      }
      for (int i_edge_node = 0; i_edge_node < 2; ++i_edge_node) {
        /* Some error checking */
        if (edge_nodes[i_edge_node].entity_dim == 2) {
          t8_global_errorf ("Error: Node %li should lie on a vertex or an edge, "
                            "but it lies on a surface.\n",
                            edge_nodes[i_edge_node].index);
          return 0;
        }
        if (edge_nodes[i_edge_node].entity_dim == 1 && edge_nodes[i_edge_node].entity_tag != edge_geometry_tag) {
          t8_global_errorf ("Error: Node %li should lie on a specific edge, "
                            "but it lies on another edge.\n",
                            edge_nodes[i_edge_node].index);
          return 0;
        }
        if (edge_nodes[i_edge_node].entity_dim == 0) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_vertex_on_edge (edge_nodes[i_edge_node].entity_tag,
                                                                            edge_geometry_tag)) {
            t8_global_errorf ("Error: Node %li should lie on a vertex which lies on an edge, "
                              "but the vertex does not lie on that edge.\n",
                              edge_nodes[i_edge_node].index);
            return 0;
          }
        }

        /* If the node lies on a vertex we retrieve its parameter on the curve */
        if (edge_nodes[i_edge_node].entity_dim == 0) {
          cad_geometry->get_cad_manager ()->t8_geom_get_parameter_of_vertex_on_edge (
            edge_nodes[i_edge_node].entity_tag, edge_geometry_tag, edge_nodes[i_edge_node].parameters.data ());
          edge_nodes[i_edge_node].entity_dim = 1;
        }
      }

      /* Abort if the edge is a line */
      if (cad_geometry->get_cad_manager ()->t8_geom_is_line (edge_geometry_tag)) {
        continue;
      }

      edge_geometries[i_tree_edges] = edge_geometry_tag;
      tree_is_linked = 1;
      parameters[0] = edge_nodes[0].parameters[0];
      parameters[1] = edge_nodes[1].parameters[0];

      /* Corrects the parameters on the edge if it is closed to prevent disorted elements. */
      if (cad_geometry->get_cad_manager ()->t8_geom_is_edge_closed (edge_geometry_tag)) {
        t8_cmesh_correct_parameters_on_closed_geometry (1, edge_geometry_tag, 2, cad_geometry, parameters);
      }

      t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (),
                              T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_tree_edges, parameters,
                              2 * sizeof (double), 0);
    }
    /* If we have found a surface we can look for the parameters. 
             * If the edge is locked for edges on surfaces we have to skip this edge */
    else if (edge_geometry_dim == 2 && edge_geometries[i_tree_edges + num_edges] >= 0) {
      /* If the node lies on a geometry with a different dimension we try to retrieve the parameters */
      for (int i_edge_node = 0; i_edge_node < 2; ++i_edge_node) {
        /* Some error checking */
        if (edge_nodes[i_edge_node].entity_dim == 2 && edge_nodes[i_edge_node].entity_tag != edge_geometry_tag) {
          t8_global_errorf ("Error: Node %li should lie on a specific face, but it lies on another face.\n",
                            edge_nodes[i_edge_node].index);
          return 0;
        }
        if (edge_nodes[i_edge_node].entity_dim == 0) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_vertex_on_face (edge_nodes[i_edge_node].entity_tag,
                                                                            edge_geometry_tag)) {
            t8_global_errorf ("Error: Node %li should lie on a vertex which lies on a face, "
                              "but the vertex does not lie on that face.\n",
                              edge_nodes[i_edge_node].index);
            return 0;
          }
        }
        if (edge_nodes[i_edge_node].entity_dim == 1) {
          if (!cad_geometry->get_cad_manager ()->t8_geom_is_edge_on_face (edge_nodes[i_edge_node].entity_tag,
                                                                          edge_geometry_tag)) {
            t8_global_errorf ("Error: Node %li should lie on an edge which lies on a face, "
                              "but the edge does not lie on that face.\n",
                              edge_nodes[i_edge_node].index);
            return 0;
          }
        }

        /* If the node lies on a vertex we retrieve its parameters on the surface */
        if (edge_nodes[i_edge_node].entity_dim == 0) {
          cad_geometry->get_cad_manager ()->t8_geom_get_parameters_of_vertex_on_face (
            edge_nodes[i_edge_node].entity_tag, edge_geometry_tag, edge_nodes[i_edge_node].parameters.data ());
          edge_nodes[i_edge_node].entity_dim = 2;
        }
        /* If the node lies on an edge we have to do the same */
        if (edge_nodes[i_edge_node].entity_dim == 1) {
          const int num_face_nodes = t8_eclass_num_vertices[eclass];
          cad_geometry->get_cad_manager ()->t8_geom_edge_parameter_to_face_parameters (
            edge_nodes[i_edge_node].entity_tag, edge_geometry_tag, num_face_nodes,
            edge_nodes[i_edge_node].parameters[0], parameters, edge_nodes[i_edge_node].parameters.data ());
          edge_nodes[i_edge_node].entity_dim = 2;
        }
      }

      /* Abort if the edge is a line */
      if (cad_geometry->get_cad_manager ()->t8_geom_is_line (edge_geometry_tag)) {
        continue;
      }

      edge_geometries[i_tree_edges + t8_eclass_num_edges[eclass]] = edge_geometry_tag;
      tree_is_linked = 1;
      parameters[0] = edge_nodes[0].parameters[0];
      parameters[1] = edge_nodes[0].parameters[1];
      parameters[2] = edge_nodes[1].parameters[0];
      parameters[3] = edge_nodes[1].parameters[1];

      /* Corrects the parameters on the surface if it is closed to prevent disorted elements. */
      for (int param_dim = 0; param_dim < 2; ++param_dim) {
        if (cad_geometry->get_cad_manager ()->t8_geom_is_surface_closed (edge_geometry_tag, param_dim)) {
          t8_cmesh_correct_parameters_on_closed_geometry (2, edge_geometry_tag, 2, cad_geometry, parameters);
        }
      }

      t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (),
                              T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_tree_edges, parameters,
                              4 * sizeof (double), 0);
    }
  }
  /* Remove the -1 used to lock the edges */
  for (int i_edge = 0; i_edge < T8_ECLASS_MAX_EDGES * 2; ++i_edge) {
    if (edge_geometries[i_edge] < 0) {
      edge_geometries[i_edge] = 0;
    }
  }
  t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, face_geometries,
                          num_faces * sizeof (int), 0);
  t8_cmesh_set_attribute (cmesh, tree_count, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edge_geometries,
                          2 * num_edges * sizeof (int), 0);

  /* Now we set the tree geometry according to the tree linkage status. */
  if (tree_is_linked) {
    t8_cmesh_set_tree_geometry (cmesh, tree_count, cad_geometry_base);
    t8_debugf ("Registering tree %li with geometry %s \n", tree_count, cad_geometry_base->t8_geom_get_name ().c_str ());
  }
  else {
    t8_cmesh_set_tree_geometry (cmesh, tree_count, linear_geometry_base);
    t8_debugf ("Registering tree %li with geometry %s \n", tree_count,
               linear_geometry_base->t8_geom_get_name ().c_str ());
  }
  return 1;
}
#endif /* T8_ENABLE_OCC */

/* fp should be set after the Nodes section, right before the tree section.
 * If vertex_indices is not NULL, it is allocated and will store
 * for each tree the indices of its vertices.
 * They are stored as arrays of long ints. 
 * If cad geometry is used, the geometry is passed as a pointer here.
 * We cannot access this geometry over the cmesh interface since the cmesh
 * is not committed yet. */
static std::optional<t8_msh_tree_vertex_indices>
t8_cmesh_msh_file_4_read_eles (t8_cmesh_t cmesh, FILE *fp, const t8_msh_node_table vertices, const int dim,
                               const t8_geometry_c *linear_geometry_base, const int use_cad_geometry,
                               [[maybe_unused]] const t8_geometry_c *cad_geometry_base, const bool store_node_data)
{
  char *line = (char *) malloc (1024), *line_modify;
  char first_word[2048] = "\0";
  size_t linen = 1024;
  t8_locidx_t tree_loop, block_loop, num_trees;
  t8_gloidx_t tree_count;
  t8_eclass_t eclass;
  t8_msh_file_node Node;
#if T8_ENABLE_OCC
  std::array<t8_msh_file_node, T8_ECLASS_MAX_CORNERS_2D> face_nodes;
  std::array<t8_msh_file_node, 2> edge_nodes;
#endif /* T8_ENABLE_OCC */
  long lnum_trees, lnum_blocks, entity_tag;
  int retval;
  int ele_type;
  int entity_dim;
  long num_ele_in_block;
  std::array<t8_msh_file_node, T8_ECLASS_MAX_CORNERS> tree_nodes;
  std::array<double, T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM> tree_vertices;
  t8_gloidx_t global_id_of_node[T8_ECLASS_MAX_CORNERS];

  T8_ASSERT (fp != NULL);
  /* Search for the line beginning with "$Elements" */
  while (!feof (fp) && strcmp (first_word, "$Elements")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2047s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num trees.\n");
      t8_debugf ("The line is %s", line);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
  }

  /* Read the line containing the number of blocks and trees */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Since t8_locidx_t could be int32 or int64, we first read the
   * number of trees in a long int and store it as t8_locidx_t later. */
  retval = sscanf (line, "%li %li %*i %*i", &lnum_blocks, &lnum_trees);
  /* Checking for read/write error */
  if (retval != 2) {
    t8_global_errorf ("Premature end of line while reading num trees and num blocks.\n");
    t8_debugf ("The line is %s", line);
    free (line);
    t8_cmesh_destroy (&cmesh);
    return std::nullopt;
  }
  num_trees = lnum_trees;
  /* Check for type conversion error */
  T8_ASSERT (num_trees == lnum_trees);

  /* Reserve memory for vertex indices */
  t8_msh_tree_vertex_indices vertex_indices (num_trees);

  tree_count = 0; /* The index of the next tree to insert */
  for (block_loop = 0; block_loop < lnum_blocks; block_loop++) {
    /* The line describing the block information looks like
     * entityDim entityTag elementType numElementsInBlock */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Premature end of line while reading trees.\n");
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
    retval = sscanf (line, "%i %li %i %li", &entity_dim, &entity_tag, &ele_type, &num_ele_in_block);
    /* Checking for read/write error */
    if (retval != 4) {
      t8_global_errorf ("Error while reading element block information.\n");
      t8_debugf ("The line is %s", line);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
    /* Check if the tree type is supported */
    if (ele_type > T8_NUM_GMSH_ELEM_CLASSES || ele_type < 0
        || t8_msh_tree_type_to_eclass[ele_type] == T8_ECLASS_COUNT) {
      t8_global_errorf ("tree type %i is not supported by t8code.\n", ele_type);
      free (line);
      t8_cmesh_destroy (&cmesh);
      return std::nullopt;
    }
    eclass = t8_msh_tree_type_to_eclass[ele_type];
    T8_ASSERT (eclass != T8_ECLASS_COUNT);

    if (t8_eclass_to_dimension[eclass] > dim) {
      t8_errorf (
        "Warning: Encountered element which dimension is greater than %d. Did you set the correct dimension?\n", dim);
    }

    /* Check if the tree is of the correct dimension */
    if (t8_eclass_to_dimension[eclass] != dim) {
      /* The trees in this block are not of the correct dimension.
       * Thus, we skip them. */
      for (tree_loop = 0; tree_loop < num_ele_in_block; tree_loop++) {
        retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
        if (retval < 0) {
          t8_global_errorf ("Premature end of line while reading trees.\n");
          free (line);
          t8_cmesh_destroy (&cmesh);
          return std::nullopt;
        }
      }
    }
    else {
      for (tree_loop = 0; tree_loop < num_ele_in_block; tree_loop++, tree_count++) {
        /* Read the next line containing tree information */
        retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
        if (retval < 0) {
          t8_global_errorf ("Premature end of line while reading trees.\n");
          free (line);
          t8_cmesh_destroy (&cmesh);
          return std::nullopt;
        }
        t8_cmesh_set_tree_class (cmesh, tree_count, eclass);

        /* The line describing the tree looks like
         * tree_number(every ele type has its own numeration) Node_1 ... Node_m
         *
         * We ignore the tree number and read the nodes.
         */
        line_modify = line;
        T8_ASSERT (strcmp (line_modify, "\0"));
        /* move line_modify to the next word in the line */
        (void) strsep (&line_modify, " ");

        /* At this point line_modify contains only the node indices. */
        const int num_nodes = t8_eclass_num_vertices[eclass];
        std::vector<t8_gloidx_t> node_indices (num_nodes, -1);
        for (int i_node = 0; i_node < num_nodes; i_node++) {
          const int t8_vertex_num = t8_msh_tree_vertex_to_t8_vertex_num[eclass][i_node];
          T8_ASSERT (strcmp (line_modify, "\0"));
          retval = sscanf (line_modify, "%li", &node_indices[t8_vertex_num]);
          if (retval != 1) {
            t8_global_errorf ("Premature end of line while reading tree.\n");
            t8_debugf ("The line is %s", line);
            free (line);
            t8_cmesh_destroy (&cmesh);
            return std::nullopt;
          }
          /* Get node from the hashtable */
          Node.index = node_indices[t8_vertex_num];
          tree_nodes[t8_vertex_num] = *vertices.find (Node);

          /* move line_modify to the next word in the line */
          (void) strsep (&line_modify, " ");
        }

        /* Add the node indices to return vector. */
        vertex_indices[tree_count] = std::move (node_indices);

        /* Set the vertices of this tree (can be removed with negative volume check) */
        for (int i_node = 0; i_node < num_nodes; i_node++) {
          tree_vertices[3 * i_node] = tree_nodes[i_node].coordinates[0];
          tree_vertices[3 * i_node + 1] = tree_nodes[i_node].coordinates[1];
          tree_vertices[3 * i_node + 2] = tree_nodes[i_node].coordinates[2];
        }
        /* Detect and correct negative volumes */
        if (t8_cmesh_tree_vertices_negative_volume (eclass, tree_vertices.data (), num_nodes)) {
          /* The volume described is negative. We need to change vertices.
           * For tets we switch 0 and 3.
           * For prisms we switch 0 and 3, 1 and 4, 2 and 5.
           * For hexahedra we switch 0 and 4, 1 and 5, 2 and 6, 3 and 7.
           * For pyramids we switch 0 and 4 */
          t8_msh_file_node temp_node;
          int num_switches = 0;
          int switch_indices[4] = { 0 };
          int iswitch;
          T8_ASSERT (t8_eclass_to_dimension[eclass] > 1);
          t8_debugf ("Correcting negative volume of tree %li\n", static_cast<long> (tree_count));
          switch (eclass) {
          case T8_ECLASS_TRIANGLE:
          case T8_ECLASS_QUAD:
            /* We switch vertex 1 and vertex 2. */
            num_switches = 2;
            switch_indices[0] = 0;
            switch_indices[1] = 2;
            break;
          case T8_ECLASS_TET:
            /* We switch vertex 0 and vertex 3. */
            num_switches = 1;
            switch_indices[0] = 3;
            break;
          case T8_ECLASS_PRISM:
            num_switches = 3;
            switch_indices[0] = 3;
            switch_indices[1] = 4;
            switch_indices[2] = 5;
            break;
          case T8_ECLASS_HEX:
            num_switches = 4;
            switch_indices[0] = 4;
            switch_indices[1] = 5;
            switch_indices[2] = 6;
            switch_indices[3] = 7;
            break;
          case T8_ECLASS_PYRAMID:
            num_switches = 1;
            switch_indices[0] = 4;
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }

          for (iswitch = 0; iswitch < num_switches; ++iswitch) {
            /* We switch node 0 + iswitch and node switch_indices[iswitch] */
            temp_node = tree_nodes[iswitch];
            tree_nodes[iswitch] = tree_nodes[switch_indices[iswitch]];
            tree_nodes[switch_indices[iswitch]] = temp_node;
          }
        }

        /* Set the vertices and global indices of this tree */
        for (int i_node = 0; i_node < num_nodes; i_node++) {
          tree_vertices[3 * i_node] = tree_nodes[i_node].coordinates[0];
          tree_vertices[3 * i_node + 1] = tree_nodes[i_node].coordinates[1];
          tree_vertices[3 * i_node + 2] = tree_nodes[i_node].coordinates[2];
          global_id_of_node[i_node] = tree_nodes[i_node].index % vertices.size ();
        }
        t8_cmesh_set_tree_vertices (cmesh, tree_count, tree_vertices.data (), num_nodes);
        t8_cmesh_set_global_vertices_of_tree (cmesh, tree_count, global_id_of_node, num_nodes);

        /* Add two arrays if store_node_data is true. One with the dimension and indices of the nodes and one with the coordinates. */
        if (store_node_data) {
          t8_store_element_node_data (cmesh, tree_count, &tree_nodes, num_nodes);
        }

        if (!use_cad_geometry) {
          /* Set the geometry of the tree to be linear.
           * If we use an cad geometry, we set the geometry in accordance,
           * if the tree is linked to a geometry or not */
          t8_cmesh_set_tree_geometry (cmesh, tree_count, linear_geometry_base);
        }
        else {
#if T8_ENABLE_OCC
          if (!t8_cmesh_process_tree_geometry (cmesh, eclass, dim, tree_count, cad_geometry_base, linear_geometry_base,
                                               tree_nodes, face_nodes, edge_nodes)) {
            free (line);
            t8_cmesh_destroy (&cmesh);
            return std::nullopt;
          }

#else  /* T8_ENABLE_OCC */
          SC_ABORTF ("OCC not linked.");
#endif /* T8_ENABLE_OCC */
        }
      }
    }
  }
  free (line);
  if (tree_count == 0) {
    t8_global_errorf ("Warning: No %iD elements found in msh file.\n", dim);
    t8_cmesh_destroy (&cmesh);
    return std::nullopt;
  }
  return std::make_optional<t8_msh_tree_vertex_indices> (vertex_indices);
}

/** This struct stores all information associated to a tree's face.
 * We need it to find neighbor trees.
 */
typedef struct
{
  t8_locidx_t ltree_id; /**< The local id of the tree this face belongs to */
  int8_t face_number;   /**< The number of that face within the tree */
  int num_vertices;     /**< The number of vertices of this face. */
  long *vertices;       /**< The indices of these vertices. */
} t8_msh_file_face_t;

/* Hash a face. The hash value is the sum of its vertex indices */
static unsigned
t8_msh_file_face_hash (const void *face, [[maybe_unused]] const void *data)
{
  t8_msh_file_face_t *Face;
  int iv;
  long sum = 0;

  Face = (t8_msh_file_face_t *) face;
  for (iv = 0; iv < Face->num_vertices; iv++) {
    sum += Face->vertices[iv];
  }
  T8_ASSERT (sum >= 0);
  return (unsigned) sum;
}

/* Two face are considered equal if they have the same vertices up
 * to renumeration. */
static int
t8_msh_file_face_equal (const void *facea, const void *faceb, [[maybe_unused]] const void *data)
{
  int iv, jv, ret;
  long vertex;
  t8_msh_file_face_t *Face_a, *Face_b;

  Face_a = (t8_msh_file_face_t *) facea;
  Face_b = (t8_msh_file_face_t *) faceb;
  /* If both have different number of vertices they can't be equal */
  ret = Face_a->num_vertices == Face_b->num_vertices;
  if (!ret) {
    return 0;
  }
  /* Check for each vertex of Face_a whether it is a vertex
   * of Face_b */
  for (iv = 0; iv < Face_a->num_vertices; iv++) {
    vertex = Face_a->vertices[iv];
    ret = 0;
    for (jv = 0; jv < Face_b->num_vertices; jv++) {
      ret |= vertex == Face_b->vertices[jv];
    }
    /* if vertex was a vertex of Face_b then ret is true here */
    if (!ret) {
      return 0;
    }
  }
  return 1;
}

/* We use this function in a loop over all elements
 * in the hash table, to free the memory of the vertices array */
static int
t8_msh_file_face_free (void **face, [[maybe_unused]] const void *data)
{
  t8_msh_file_face_t *Face;

  Face = *(t8_msh_file_face_t **) face;
  T8_FREE (Face->vertices);
  return 1;
}

/* Given a face and a cmesh set the face as a domain boundary.
 * We use this function in a loop over all hashed faces.
 */
static int
t8_msh_file_face_set_boundary (void **face, const void *data)
{
  t8_msh_file_face_t *Face;
  t8_cmesh_t cmesh;
  t8_gloidx_t gtree_id;

  cmesh = (t8_cmesh_t) data;
  Face = *(t8_msh_file_face_t **) face;

  /* Get the global tree id */
  gtree_id = Face->ltree_id;
  /* Set the Face as a domain boundary */
  t8_cmesh_set_join (cmesh, gtree_id, gtree_id, Face->face_number, Face->face_number, 0);
  return 1;
}

/* Given two faces and the classes of their volume trees,
 * compute the orientation of the faces to each other */
static int
t8_msh_file_face_orientation (const t8_msh_file_face_t *Face_a, const t8_msh_file_face_t *Face_b,
                              const t8_eclass_t tree_class_a, const t8_eclass_t tree_class_b)
{
  long vertex_zero; /* The number of the first vertex of the smaller face */
  const t8_msh_file_face_t *smaller_Face, *bigger_Face;
  int compare, iv;
  t8_eclass_t bigger_class;
  int orientation = -1;

  compare = t8_eclass_compare (tree_class_a, tree_class_b);
  if (compare > 0) {
    /* tree_class_a is bigger than tree_class_b */
    smaller_Face = Face_b;
    bigger_Face = Face_a;
    bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
  }
  else if (compare < 0) {
    /* tree_class_a is smaller than tree_class_b */
    smaller_Face = Face_a;
    bigger_Face = Face_b;
    bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
  }
  else {
    /* both classes are the same, thus
     * the face with the smaller face id is the smaller one */
    if (Face_a->face_number < Face_b->face_number) {
      smaller_Face = Face_a;
      bigger_Face = Face_b;
      bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
    }
    else {
      smaller_Face = Face_b;
      bigger_Face = Face_a;
      bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
    }
  }
  vertex_zero = smaller_Face->vertices[0];
  /* Find which point in the bigger face is vertex_zero */
  for (iv = 0; iv < t8_eclass_num_vertices[bigger_class]; iv++) {
    if (vertex_zero == bigger_Face->vertices[iv]) {
      /* We found the corresponding vertex */
      orientation = iv;
      /* set condition to break the loop */
      iv = t8_eclass_num_vertices[bigger_class];
    }
  }
  T8_ASSERT (orientation >= 0);
  return orientation;
}

/* Given the number of vertices and for each element a list of its
 * vertices, find the neighborship relations of each element */
/* This routine does only find neighbors between local trees.
 * Use with care if cmesh is partitioned. */
static void
t8_cmesh_msh_file_find_neighbors (t8_cmesh_t cmesh, const t8_msh_tree_vertex_indices vertex_indices)
{
  sc_hash_t *faces;
  t8_msh_file_face_t *Face, **pNeighbor, *Neighbor;
  sc_mempool_t *face_mempool;
  t8_gloidx_t gtree_it;
  t8_gloidx_t gtree_id, gtree_neighbor;
  t8_eclass_t eclass, face_class, neighbor_tclass;
  int num_face_vertices, face_it, vertex_it;
  int retval, orientation;
  t8_stash_class_struct_t *class_entry;

  face_mempool = sc_mempool_new (sizeof (t8_msh_file_face_t));
  faces = sc_hash_new (t8_msh_file_face_hash, t8_msh_file_face_equal, cmesh, NULL);

  /* TODO: Does currently not work with partitioned cmesh */
  T8_ASSERT (!cmesh->set_partition);
  /* The cmesh is not allowed to be committed yet */
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  t8_debugf ("Starting to find tree neighbors\n");
  /* Iterate over all local trees */
  for (gtree_it = 0; gtree_it < (t8_gloidx_t) cmesh->stash->classes.elem_count; gtree_it++) {
    /* We get the class of the current tree.
     * Since we know that the trees were put into the stash in order
     * of their tree id's, we can just read the corresponding entry from
     * the stash.
     * !WARNING: This does not work in general to find the class of a tree
     *    since the order in which the trees are added to the stash is arbitrary.
     *    Use t8_stash_class_bsearch in tat case.
     */
    class_entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx (&cmesh->stash->classes, gtree_it);
    T8_ASSERT (class_entry->id == gtree_it);
    eclass = class_entry->eclass;
    /* Get the vertices of that tree */
    auto tree_vertices = vertex_indices[gtree_it];
    /* loop over all faces of the tree */
    for (face_it = 0; face_it < t8_eclass_num_faces[eclass]; face_it++) {
      /* Create the Face struct */
      Face = (t8_msh_file_face_t *) sc_mempool_alloc (face_mempool);
      /* Get its eclass and the number of vertices */
      face_class = (t8_eclass_t) t8_eclass_face_types[eclass][face_it];
      num_face_vertices = t8_eclass_num_vertices[face_class];
      Face->vertices = T8_ALLOC (long, num_face_vertices);
      Face->num_vertices = num_face_vertices;
      Face->ltree_id = gtree_it;
      Face->face_number = face_it;
      /* Copy the vertices of the face to the face struct */
      for (vertex_it = 0; vertex_it < num_face_vertices; vertex_it++) {
        Face->vertices[vertex_it] = tree_vertices[t8_face_vertex_to_tree_vertex[eclass][face_it][vertex_it]];
      }
      /* Try to insert the face into the hash */
      retval = sc_hash_insert_unique (faces, Face, (void ***) &pNeighbor);
      if (!retval) {
        /* The face was already in the hash */
        Neighbor = *pNeighbor;
        T8_ASSERT (Neighbor->ltree_id != gtree_it);
        /* The current tree is a neighbor to the tree Neighbor->ltree_id */
        /* We need to identify the face number and the orientation */
        gtree_id = gtree_it;
        gtree_neighbor = Neighbor->ltree_id;
        /* Compute the orientation of the face connection */
        /* Get the element class of the neighbor tree */
        class_entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx (&cmesh->stash->classes, Neighbor->ltree_id);
        T8_ASSERT (class_entry->id == Neighbor->ltree_id);
        neighbor_tclass = class_entry->eclass;
        /* Calculate the orientation */
        orientation = t8_msh_file_face_orientation (Face, Neighbor, eclass, neighbor_tclass);
        /* Set the face connection */
        t8_cmesh_set_join (cmesh, gtree_id, gtree_neighbor, face_it, Neighbor->face_number, orientation);
        T8_FREE (Face->vertices);
        sc_mempool_free (face_mempool, Face);
        /* We can free the neighbor here as well */
        sc_hash_remove (faces, Neighbor, NULL);
        T8_FREE (Neighbor->vertices);
        sc_mempool_free (face_mempool, Neighbor);
      }
    }
  }
  /* The remaining faces are domain boundaries */
  sc_hash_foreach (faces, t8_msh_file_face_set_boundary);
  /* Free the faces and the hash */
  sc_hash_foreach (faces, t8_msh_file_face_free);
  sc_mempool_destroy (face_mempool);
  sc_hash_destroy (faces);
  t8_debugf ("Done finding tree neighbors.\n");
}

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

/* This is a helper function to properly register the 
 * geometries for the cmesh created in t8_cmesh_from_msh_file.
 * It should be called by all processes of the cmesh.
 * Returns 1 on success, 0 on cad usage error: use_cad_geometry true, but OCC not linked.
 * The linear_geometry pointer will point to the newly created linear geometry.
 * The cad_geometry pointer will point to the newly created cad geometry, or to NULL if
 * no cad geometry is used.
 */
static int
t8_cmesh_from_msh_file_register_geometries (t8_cmesh_t cmesh, const int use_cad_geometry,
                                            [[maybe_unused]] const char *fileprefix,
                                            const t8_geometry_c **linear_geometry, const t8_geometry_c **cad_geometry)
{
  /* Register linear geometry */
  *linear_geometry = t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  if (use_cad_geometry) {
#if T8_ENABLE_OCC
    *cad_geometry = t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, std::string (fileprefix));
#else /* !T8_ENABLE_OCC */
    *cad_geometry = NULL;
    return 0;
#endif
  }
  return 1;
}

t8_cmesh_t
t8_cmesh_from_msh_file (const char *fileprefix, const int partition, sc_MPI_Comm comm, const int dim,
                        const int main_proc, const int use_cad_geometry)
{
  int mpirank, mpisize, mpiret;
  t8_cmesh_t cmesh;
  char current_file[BUFSIZ];
  FILE *file;
  t8_gloidx_t num_trees, first_tree, last_tree = -1;
  int main_proc_read_successful = 0;
  int msh_version;
  const t8_geometry_c *cad_geometry = NULL;
  const t8_geometry_c *linear_geometry = NULL;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* TODO: implement partitioned input using gmesh's
   * partitioned files.
   * Or using a single file and computing the partition on the run. */
  T8_ASSERT (partition == 0 || (main_proc >= 0 && main_proc < mpisize));

  /* initialize cmesh structure */
  t8_cmesh_init (&cmesh);
  /* Setting the dimension by hand is necessary for partitioned
   * commit, since there are process without any trees. So the cmesh would
   * not know its dimension on these processes. */
  t8_cmesh_set_dimension (cmesh, dim);

  /* Register the geometries for the cmesh. */
  const int registered_geom_success
    = t8_cmesh_from_msh_file_register_geometries (cmesh, use_cad_geometry, fileprefix, &linear_geometry, &cad_geometry);
  if (!registered_geom_success) {
    /* Registering failed */
    t8_errorf ("OCC is not linked. Cannot use cad geometry.\n");
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }

  if (!partition || mpirank == main_proc) {
    snprintf (current_file, BUFSIZ, "%s.msh", fileprefix);
    /* Open the file */
    t8_debugf ("Opening file %s\n", current_file);
    file = fopen (current_file, "r");
    if (file == NULL) {
      t8_global_errorf ("Could not open file %s\n", current_file);
      t8_cmesh_destroy (&cmesh);

      if (partition) {
        /* Communicate to the other processes that reading failed. */
        main_proc_read_successful = 0;
        sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
      }
      return NULL;
    }
    /* Check if msh-file version is compatible. */
    msh_version = t8_cmesh_check_version_of_msh_file (file);
    if (msh_version < 1) {
      /* If reading the MeshFormat-number failed or the version is incompatible, close the file */
      fclose (file);
      t8_debugf ("The reading process of the msh-file has failed and the file has been closed.\n");
      t8_cmesh_destroy (&cmesh);

      if (partition) {
        /* Communicate to the other processes that reading failed. */
        main_proc_read_successful = 0;
        sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
      }
      return NULL;
    }
    /* read nodes from the file */
    std::optional<t8_msh_tree_vertex_indices> indices;
    switch (msh_version) {
    case 2: {
      if (use_cad_geometry) {
        fclose (file);
        t8_errorf ("WARNING: The cad geometry is only supported for msh files of version 4\n");
        t8_cmesh_destroy (&cmesh);
        if (partition) {
          /* Communicate to the other processes that reading failed. */
          main_proc_read_successful = 0;
          sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
        }
        return NULL;
      }
      auto vertices_opt = t8_msh_file_2_read_nodes (file);
      if (!vertices_opt) {
        fclose (file);
        t8_cmesh_destroy (&cmesh);
        if (partition) {
          /* Communicate to the other processes that reading failed. */
          main_proc_read_successful = 0;
          sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
        }
        return NULL;
      }
      indices = t8_cmesh_msh_file_2_read_eles (cmesh, file, *vertices_opt, dim);
      break;
    }

    case 4: {
      auto vertices_opt = t8_msh_file_4_read_nodes (file);
      if (!vertices_opt) {
        fclose (file);
        t8_cmesh_destroy (&cmesh);
        if (partition) {
          /* Communicate to the other processes that reading failed. */
          main_proc_read_successful = 0;
          sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
        }
        return NULL;
      }
      indices = t8_cmesh_msh_file_4_read_eles (cmesh, file, *vertices_opt, dim, linear_geometry, use_cad_geometry,
                                               cad_geometry, true);
      break;
    }

    default:
      break;
    }
    /* close the file and free the memory for the nodes */
    fclose (file);
    if (!indices) {
      t8_cmesh_destroy (&cmesh);
      if (partition) {
        /* Communicate to the other processes that reading failed. */
        main_proc_read_successful = 0;
        sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
      }
      return NULL;
    }
    else
      t8_cmesh_msh_file_find_neighbors (cmesh, *indices);

    main_proc_read_successful = 1;
  }

  if (partition) {
    /* Communicate whether main proc read the cmesh successful.
     * If the main process failed then it called this Bcast already and
     * terminated. If it was successful, it calls the Bcast now. */
    sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
    if (!main_proc_read_successful) {
      t8_debugf ("Main process could not read cmesh successfully.\n");
      t8_cmesh_destroy (&cmesh);
      return NULL;
    }
    /* The cmesh is not yet committed, since we set the partitioning before */
    if (mpirank == main_proc) {
      /* The main_proc process sends the number of trees to
       * all processes. This is used to fill the partition table
       * that says that all trees are on main_proc and zero on everybody else. */
      num_trees = cmesh->stash->classes.elem_count;
      first_tree = 0;
      last_tree = num_trees - 1;
      T8_ASSERT (cmesh->dimension == dim);
    }
    /* bcast the global number of trees */
    sc_MPI_Bcast (&num_trees, 1, T8_MPI_GLOIDX, main_proc, comm);
    /* Set the first and last trees on this rank.
     * No rank has any trees except the main_proc */
    if (mpirank < main_proc) {
      first_tree = 0;
      last_tree = -1;
    }
    else if (mpirank > main_proc) {
      first_tree = num_trees;
      last_tree = num_trees - 1;
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }

  /* Commit the cmesh */
  T8_ASSERT (cmesh != NULL);
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  return cmesh;
}

T8_EXTERN_C_END ();
