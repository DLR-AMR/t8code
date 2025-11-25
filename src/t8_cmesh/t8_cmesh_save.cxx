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

/** \file t8_cmesh_save.cxx
 *
 * We define routines to save and load a cmesh to/from the file system.
 *
 * TODO: document this file
 */

#include <t8_version.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_save.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_base.h>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

/**
 * This macro is called to check a condition and if not fulfilled
 * close the file and exit the function
 */
#define T8_SAVE_CHECK_CLOSE(x, fp) \
  if (!(x)) { \
    t8_errorf ("file i/o error. Condition %s not fulfilled. " \
               "Line %i in file %s\n", \
               #x, __LINE__, __FILE__); \
    fclose (fp); \
    return 0; \
  }

/* Write the neighbor data of all ghosts */
static int
t8_cmesh_save_ghost_neighbors (const t8_cmesh_t cmesh, FILE *fp)
{
  t8_locidx_t ighost;
  t8_cghost_t ghost;
  t8_gloidx_t *face_neigh;
  int8_t *ttf;
  int iface, num_faces;
  int ret;
  char buffer[BUFSIZ];

  ret = fprintf (fp, "\n--- Ghost neighbor section ---\n\n");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    /* Reset the buffer */
    buffer[0] = '\0';
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees, ighost, &face_neigh, &ttf);
    /* Iterate over all faces to write the face neighbor information */
    num_faces = t8_eclass_num_faces[ghost->eclass];
    for (iface = 0; iface < num_faces; iface++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%lli %i%s", (long long) face_neigh[iface],
                ttf[iface], iface == num_faces - 1 ? "" : ", ");
    }
    ret = fprintf (fp, "%s\n", buffer);
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  }
  return 1;
}

static int
t8_cmesh_load_ghost_attributes (const t8_cmesh_t cmesh, FILE *fp)
{

  t8_locidx_t ighost;
  t8_cghost_t ghost;
  t8_gloidx_t *face_neigh;
  int8_t *ttf;
  int iface, num_faces, ttf_entry;
  long long neigh;
  int ret;

  ret = fscanf (fp, "\n--- Ghost neighbor section ---\n\n");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);

  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    /* Get a pointer to the ghost */
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees, ighost, &face_neigh, &ttf);
    /* Read the neighbor information */
    num_faces = t8_eclass_num_faces[ghost->eclass];
    for (iface = 0; iface < num_faces; iface++) {
      ret = fscanf (fp, "%lli %i%*c", &neigh, &ttf_entry);
      T8_SAVE_CHECK_CLOSE (ret == 2, fp);
      ttf[iface] = ttf_entry;
      face_neigh[iface] = neigh;
    }
  }
  return 1;
}

static int
t8_cmesh_save_ghosts (const t8_cmesh_t cmesh, FILE *fp)
{
  t8_locidx_t ighost;
  t8_cghost_t ghost;
  int ret;

  ret = fprintf (fp, "\n--- Ghost section ---");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  ret = fprintf (fp, "\n\nGhosts %li\n\n", (long) cmesh->num_ghosts);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    ghost = t8_cmesh_trees_get_ghost (cmesh->trees, ighost);
    ret = fprintf (fp, "treeid %lli\neclass %i\n\n", (long long) ghost->treeid, (int) ghost->eclass);
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  }
  return 1;
}

static int
t8_cmesh_load_ghosts (const t8_cmesh_t cmesh, FILE *fp)
{
  t8_locidx_t ighost;
  t8_cghost_t ghost;
  int ret;
  int eclass;
  long num_ghosts;
  long long global_id;

  ret = fscanf (fp, "\n--- Ghost section ---\n");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  /* Check whether the number of ghosts is correct */
  ret = fscanf (fp, "Ghosts %li\n", &num_ghosts);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  T8_SAVE_CHECK_CLOSE (num_ghosts == cmesh->num_ghosts, fp);
  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    /* Get a pointer to the ghost */
    ghost = t8_cmesh_trees_get_ghost (cmesh->trees, ighost);
    /* Read the ghost's global id */
    ret = fscanf (fp, "treeid %lli\n", &global_id);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    T8_SAVE_CHECK_CLOSE (0 <= global_id && global_id < cmesh->num_trees, fp);
    ghost->treeid = global_id;
    /* read the ghost's eclass */
    ret = fscanf (fp, "eclass %i\n", &eclass);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    ghost->eclass = (t8_eclass_t) eclass;
    /* Ignore the neighbor offset */
    ret = fscanf (fp, "neigh_offset %*i\n");
    T8_SAVE_CHECK_CLOSE (ret == 0, fp);
    ret = fscanf (fp, "\n");
    T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  }
  return 1;
}

/* Load all attributes that were stored in a file */
static int
t8_cmesh_load_tree_attributes (const t8_cmesh_t cmesh, FILE *fp)
{
  double *vertices = NULL;
  t8_locidx_t itree;
  long treeid, neighbor;
  t8_ctree_t tree;
  int att, i, num_vertices, num_faces, iface;
  int ret, ttf_entry;
  t8_locidx_t *face_neighbors;
  int8_t *ttf;
  t8_stash_attribute_struct_t att_struct;

  ret = fscanf (fp, "\n--- Tree attribute section ---\n");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  /* loop over all trees */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbors, &ttf);
    /* Allocate memory for the temporary vertices array.
     * This memory is reused for each tree */
    num_vertices = t8_eclass_num_vertices[tree->eclass];
    ret = fscanf (fp, "tree %li\n", &treeid);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    /* Load the tree neighbors */
    num_faces = t8_eclass_num_faces[tree->eclass];
    ret = fscanf (fp, "Neighbors:");
    T8_SAVE_CHECK_CLOSE (ret == 0, fp);
    for (iface = 0; iface < num_faces; iface++) {
      /* Read the face neighbors tree id and the tree_to_face entry */
      ret = fscanf (fp, "%li %i%*c", &neighbor, &ttf_entry);
      T8_SAVE_CHECK_CLOSE (ret == 2, fp);
      face_neighbors[iface] = neighbor;
      ttf[iface] = ttf_entry;
    }
    vertices = SC_REALLOC (vertices, double, 3 * num_vertices);
    for (att = 0; att < tree->num_attributes; att++) {
      /* Loop over all attributes of this tree that we need to read */
      /* Check whether this attribute really belongs to tree */
      T8_SAVE_CHECK_CLOSE (treeid == (long) itree, fp);
      ret = fscanf (fp, "id %i\nkey %i\n", &att_struct.package_id, &att_struct.key);
      T8_SAVE_CHECK_CLOSE (ret == 2, fp);
      /* We currently only support vertices as attributes.
       * Those have t8 package id and key 0 */
      T8_SAVE_CHECK_CLOSE (att_struct.package_id > 0 && att_struct.key == 0, fp);
      /* TODO: We set the package id to match the one of t8code manually.
       * See the comment above. */
      att_struct.package_id = t8_get_package_id ();

      /* read the size of the attribute */
      ret = fscanf (fp, "size %zu\n", &att_struct.attr_size);
      T8_SAVE_CHECK_CLOSE (ret == 1, fp);
      /* Read the vertices */
      for (i = 0; i < num_vertices; i++) {
        ret = fscanf (fp, "%lf %lf %lf\n", vertices + 3 * i, vertices + 3 * i + 1, vertices + 3 * i + 2);
        T8_SAVE_CHECK_CLOSE (ret == 3, fp);
      }
      att_struct.attr_data = vertices;
      att_struct.is_owned = 0;
      att_struct.id = itree + cmesh->first_tree;
      /* Now we read the attribute and can add it to the tree */
      t8_cmesh_trees_add_attribute (cmesh->trees, 0, &att_struct, itree, 0);
    }
  }
  SC_FREE (vertices);
  return 1;
}

static int
t8_cmesh_save_tree_attribute (const t8_cmesh_t cmesh, FILE *fp)
{
  double *vertices;
  t8_locidx_t itree;
  t8_ctree_t tree;
  int num_vertices;
  int ret, i;
  size_t att_size;
  t8_locidx_t *face_neigh;
  int8_t *ttf;
  int num_faces;
  char buffer[BUFSIZ] = "";

  ret = fprintf (fp, "\n--- Tree attribute section ---\n");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  /* For each tree, write its attribute */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    /* TODO: Currently we only support storing of the tree vertices as only attribute.
     *        This attribute has package id t8_get_package_id() and key 0.
     */
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neigh, &ttf);
    ret = fprintf (fp, "\ntree %li\n", (long) itree);
    /* Iterate over all faces to write the face neighbor information */
    num_faces = t8_eclass_num_faces[tree->eclass];
    for (i = 0; i < num_faces; i++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%li %i%s", (long) face_neigh[i], ttf[i],
                i == num_faces - 1 ? "" : ", ");
    }
    ret = fprintf (fp, "Neighbors: %s\n", buffer);
    /* Clear the buffer such that strlen returns 0 */
    buffer[0] = '\0';
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
    num_vertices = t8_eclass_num_vertices[tree->eclass];
    /* Write the attributes that are vertices */
    vertices = (double *) t8_cmesh_trees_get_attribute (cmesh->trees, itree, t8_get_package_id (), 0, &att_size, 0);
    if (vertices != NULL) {
      /* We have an attribute that is stored with key 0, we treat it as tree vertices */
      num_vertices = t8_eclass_num_vertices[tree->eclass];
      /* additional check that the attribute has as many bytes as the
       * coordinates for this tree */
      if (att_size != num_vertices * 3 * sizeof (double)) {
        /* TODO: We currently do not support saving of attributes different
         * to the tree vertices */
        fclose (fp);
        t8_errorf ("We do not support saving cmeshes with trees that "
                   "have attributes different to the tree vertices.\n");
        return 0;
      }
      ret = fprintf (fp, "id %i\nkey %i\n", t8_get_package_id (), 0);
      T8_SAVE_CHECK_CLOSE (ret > 0, fp);
      ret = fprintf (fp, "size %zd\n", num_vertices * 3 * sizeof (double));
      T8_ASSERT (strlen (buffer) == 0);
      for (i = 0; i < num_vertices; i++) {
        /* For each vertex, we write its three coordinates in a line */
        snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%e %e %e\n", vertices[3 * i],
                  vertices[3 * i + 1], vertices[3 * i + 2]);
      }
      fprintf (fp, "%s", buffer);
      /* Clear the buffer such that strlen returns 0 */
      buffer[0] = '\0';
    }
  }
  return 1;
}

static int
t8_cmesh_save_trees (const t8_cmesh_t cmesh, FILE *fp)
{
  t8_locidx_t itree;
  t8_ctree_t tree;
  int ret;

  T8_ASSERT (fp != NULL);
  ret = fprintf (fp, "Total bytes for trees and ghosts %zd\n", t8_cmesh_trees_size (cmesh->trees));
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  ret = fprintf (fp, "\n--- Tree section ---");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  ret = fprintf (fp, "\n\nTrees %li\n\n", (long) cmesh->num_local_trees);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* For each tree, write its metadata */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    tree = t8_cmesh_trees_get_tree (cmesh->trees, itree);
    ret = fprintf (fp, "eclass %i\n", (int) tree->eclass);
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
    if (tree->num_attributes != 1) {
      /* TODO: We currently do not support saving of attributes different
       * to the tree vertices */
      fclose (fp);
      t8_errorf ("We do not support saving cmeshes with trees that have attributes different to the tree vertices.\n");
      return 0;
    }
    ret = fprintf (fp, "num_attributes %i\nSize of attributes %zd\n\n", tree->num_attributes,
                   t8_cmesh_trees_attribute_size (tree));
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  }
  return 1;
}

/* Load all tree data (eclasses, neighbors, vertex coordinates,...)
 * from a cmesh file into a cmesh */
static int
t8_cmesh_load_trees (const t8_cmesh_t cmesh, FILE *fp)
{
  size_t bytes_for_trees, att_bytes;
  t8_locidx_t itree;
  int eclass, num_atts;
  int ret;
  long num_trees;

  ret = fscanf (fp, "Total bytes for trees and ghosts %zu\n", &bytes_for_trees);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp); /* The bytes_for_trees data is currently not used */
  ret = fscanf (fp, "\n--- Tree section ---\n");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  t8_cmesh_trees_init (&cmesh->trees, 1, cmesh->num_local_trees, cmesh->num_ghosts);
  t8_cmesh_trees_start_part (cmesh->trees, 0, 0, cmesh->num_local_trees, 0, cmesh->num_ghosts, 1);
  /* Read the number of trees to come and check for consistency */
  ret = fscanf (fp, "Trees %li\n", &num_trees);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  T8_SAVE_CHECK_CLOSE (num_trees == cmesh->num_local_trees, fp);
  /* Read all the trees from the file and add them to the
   * trees structure */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    ret = fscanf (fp, "eclass %i\n", &eclass);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    t8_cmesh_trees_add_tree (cmesh->trees, itree, 0, (t8_eclass_t) eclass);

    /* After adding the tree, we set its face neighbors and face orientation */
    (void) t8_cmesh_trees_get_tree (cmesh->trees, itree);
    /* Check whether the number of attribute is really 1 */
    ret = fscanf (fp, "num_attributes %i\nSize of attributes %zu\n", &num_atts, &att_bytes);
    T8_SAVE_CHECK_CLOSE (ret == 2, fp);
    T8_SAVE_CHECK_CLOSE (num_atts == 1, fp);
    T8_SAVE_CHECK_CLOSE (att_bytes == 3 * sizeof (double) * t8_eclass_num_vertices[eclass], fp);
    /* If there really is one attribute it must be the trees vertices and thus
     * we initialize the tree's attribute accordingly */
    t8_cmesh_trees_init_attributes (cmesh->trees, itree, 1, att_bytes);
  }
  return 1;
}

static int
t8_cmesh_save_header (const t8_cmesh_t cmesh, FILE *fp)
{
  int ret;
  int eclass;

  T8_ASSERT (fp != NULL);
  ret = fprintf (fp, "This is %s, file format version %u.\n\n", t8_get_package_string (), T8_CMESH_FORMAT);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* Write 0 for replicated and 1 for partitioned cmesh */
  ret = fprintf (fp, "Partitioned %i\n", cmesh->set_partition != 0);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* Write the rank of this process and the total number pf processes */
  ret = fprintf (fp, "Rank %i of %i\n", cmesh->mpirank, cmesh->mpisize);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* Write the dimension of the cmesh */
  ret = fprintf (fp, "dim %i\n", cmesh->dimension);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* Write the number of global and local trees and the number of ghosts */
  ret = fprintf (fp, "num_trees %lli\n", (long long) cmesh->num_trees);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  ret = fprintf (fp, "num_local_trees %li\n", (long) cmesh->num_local_trees);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  ret = fprintf (fp, "num_ghosts %li\n", (long) cmesh->num_ghosts);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);

  /* Write the number of local trees for each eclass */
  ret = fprintf (fp, "num_local_trees_per_eclass ");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; eclass++) {
    ret = fprintf (fp, "%li%s", (long) cmesh->num_local_trees_per_eclass[eclass],
                   eclass == T8_ECLASS_COUNT - 1 ? "\n" : ", ");
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  }
  /* Write the number of trees for each eclass */
  ret = fprintf (fp, "num_trees_per_eclass ");
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; eclass++) {
    ret = fprintf (fp, "%lli%s", (long long) cmesh->num_trees_per_eclass[eclass],
                   eclass == T8_ECLASS_COUNT - 1 ? "\n" : ", ");
    T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  }

  /* Write the global id of the first local tree and if the first tree is shared */
  ret = fprintf (fp, "first_tree %lli\n", (long long) cmesh->first_tree);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  ret = fprintf (fp, "first_tree_shared %i\n", cmesh->first_tree_shared);
  T8_SAVE_CHECK_CLOSE (ret > 0, fp);
  return 1;
}

/* Read the number of trees, dimension, etc. from a saved cmesh file.
 * If anything goes wrong, the file is closed and 0 is returned */
static int
t8_cmesh_load_header (const t8_cmesh_t cmesh, FILE *fp)
{
  int file_format, save_rank, save_mpisize;
  int ieclass;
  int ret;
  long long global_num_trees, tree_per_class, first_local_tree;
  long local_num_trees, local_num_ghosts;
  int first_shared;

  /* Check whether the file was saved with the current file format */
  ret = fscanf (fp, "This is %*[^ ] %*[^ ] file format version %i.\n", &file_format);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  if (file_format != T8_CMESH_FORMAT) {
    /* The file was saved with an old format and we cannot read it any more */
    t8_errorf ("Input file is in an old format that we cannot read anymore.\n");
    fclose (fp);
    return 0;
  }

  /* Read whether the cmesh is partitioned, its dimension and the number of
   * trees. Also read which MPI rank wrote the file */
  ret = fscanf (fp, "Partitioned %i\n", &cmesh->set_partition);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  ret = fscanf (fp, "Rank %i of %i\n", &save_rank, &save_mpisize);
  T8_SAVE_CHECK_CLOSE (ret == 2, fp);
  /* Check if the rank and mpisize stored were valid */
  T8_SAVE_CHECK_CLOSE (0 <= save_rank && save_rank < save_mpisize, fp);
  /* It does not make sense to load a cmesh on a rank smaller than the one that
   * saved it. */
  T8_SAVE_CHECK_CLOSE (cmesh->mpirank <= save_rank && cmesh->mpisize <= save_mpisize, fp);
  ret = fscanf (fp, "dim %i\n", &cmesh->dimension);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  /* Check if the read dimension is in the correct range */
  T8_SAVE_CHECK_CLOSE (cmesh->dimension >= 0 && cmesh->dimension <= 3, fp);
  /* Since t8_gloidx_t and t8_locidx_t are integer datatypes that are not fixed,
   * we first read the tree numbers into fixed datatypes (long long for gloidx and
   * long for locidx) */
  ret = fscanf (fp, "num_trees %lli\nnum_local_trees %li\n", &global_num_trees, &local_num_trees);
  T8_SAVE_CHECK_CLOSE (ret == 2, fp);
  T8_SAVE_CHECK_CLOSE (local_num_trees <= global_num_trees, fp);
  cmesh->num_trees = (t8_gloidx_t) global_num_trees;
  cmesh->num_local_trees = (t8_locidx_t) local_num_trees;
  /* Read the number of ghost trees */
  ret = fscanf (fp, "num_ghosts %li\n", &local_num_ghosts);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  T8_SAVE_CHECK_CLOSE (0 <= local_num_ghosts && local_num_ghosts < cmesh->num_trees, fp);
  cmesh->num_ghosts = local_num_ghosts;

  /* Read the number of local trees per eclass */
  ret = fscanf (fp, "num_local_trees_per_eclass");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  for (ieclass = T8_ECLASS_ZERO; ieclass < T8_ECLASS_COUNT; ieclass++) {
    ret = fscanf (fp, "%lli,", &tree_per_class);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    cmesh->num_local_trees_per_eclass[ieclass] = (t8_locidx_t) tree_per_class;
  }
  /* Read the number of trees per eclass */
  ret = fscanf (fp, "\nnum_trees_per_eclass");
  T8_SAVE_CHECK_CLOSE (ret == 0, fp);
  for (ieclass = T8_ECLASS_ZERO; ieclass < T8_ECLASS_COUNT; ieclass++) {
    ret = fscanf (fp, "%lli,", &tree_per_class);
    T8_SAVE_CHECK_CLOSE (ret == 1, fp);
    cmesh->num_trees_per_eclass[ieclass] = (t8_gloidx_t) tree_per_class;
  }
  /* Read the first local tree id and whether the first tree is shared */
  ret = fscanf (fp, "\nfirst_tree %lli\n", &first_local_tree);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  T8_SAVE_CHECK_CLOSE (0 <= first_local_tree && first_local_tree < global_num_trees, fp);
  cmesh->first_tree = first_local_tree;
  ret = fscanf (fp, "first_tree_shared %i\n", &first_shared);
  T8_SAVE_CHECK_CLOSE (ret == 1, fp);
  cmesh->first_tree_shared = first_shared;
  T8_SAVE_CHECK_CLOSE (!cmesh->set_partition || cmesh->first_tree_shared == 0 || cmesh->first_tree_shared == 1, fp);
  return 1;
}

int
t8_cmesh_save (const t8_cmesh_t cmesh, const char *fileprefix)
{
  FILE *fp;
  char filename[BUFSIZ];
  int has_linear_geom = 0;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (!cmesh->set_partition && cmesh->mpirank != 0) {
    /* If the cmesh is replicated, only rank 0 writes it */
    return 1;
  }

  /* Check that the only registered geometry is the linear geometry and
   * that this geometry is used for all trees. */
  if (cmesh->geometry_handler->get_num_geometries () == 1) {
    /* Get the stored geometry and the linear geometry and compare their names. */
    const t8_geometry *geom = cmesh->geometry_handler->get_unique_geometry ();
    has_linear_geom = geom->t8_geom_get_type () == T8_GEOMETRY_TYPE_LINEAR;
  }
  if (!has_linear_geom) {
    /* This cmesh does not have the linear geometry for all trees. */
    t8_errorf ("Error when saving cmesh. Cmesh has more than one geometry or the geometry is not linear.\n");
    return 0;
  }

  /* Create the output filename as fileprefix_RANK.cmesh,
   * where we write RANK with 4 significant figures. */
  snprintf (filename, BUFSIZ, "%s_%04i.cmesh", fileprefix, cmesh->mpirank);

  /* Open the file in write mode */
  fp = fopen (filename, "w");
  if (fp == NULL) {
    /* Could not open file */
    t8_errorf ("Error when opening file %s.\n", filename);
    return 0;
  }
  /* Write all metadata of the cmesh */
  if (!t8_cmesh_save_header (cmesh, fp)) {
    t8_errorf ("Error when opening file %s.\n", filename);
    return 0;
  }
  /* Write all metadata of the trees */
  if (!t8_cmesh_save_trees (cmesh, fp)) {
    t8_errorf ("Error when opening file %s.\n", filename);
    return 0;
  }
  if (cmesh->set_partition) {
    /* Write all ghost metadata */
    if (!t8_cmesh_save_ghosts (cmesh, fp)) {
      t8_errorf ("Error when opening file %s.\n", filename);
      return 0;
    }
  }
  if (!t8_cmesh_save_tree_attribute (cmesh, fp)) {
    /* Write all tree attributes */
    t8_errorf ("Error when opening file %s.\n", filename);
    return 0;
  }
  if (cmesh->set_partition) {
    /* Write all ghost neighbor data */
    if (!t8_cmesh_save_ghost_neighbors (cmesh, fp)) {
      t8_errorf ("Error when opening file %s.\n", filename);
      return 0;
    }
  }
  /* Close the file */
  fclose (fp);
  return 1;
}

#undef T8_SAVE_CHECK_CLOSE

t8_cmesh_t
t8_cmesh_load (const char *filename, sc_MPI_Comm comm)
{
  FILE *fp;
  t8_cmesh_t cmesh;
  int mpiret;

  /* Open the file in read mode */
  fp = fopen (filename, "r");
  if (fp == NULL) {
    /* Could not open file */
    t8_errorf ("Error when opening file %s.\n", filename);
    return NULL;
  }
  t8_cmesh_init (&cmesh);
  /* Read all metadata of the cmesh */
  if (!t8_cmesh_load_header (cmesh, fp)) {
    t8_errorf ("Error when opening file %s.\n", filename);
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }
  /* Read all metadata of the trees */
  if (!t8_cmesh_load_trees (cmesh, fp)) {
    t8_errorf ("Error when opening file %s.\n", filename);
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }
  if (cmesh->set_partition) {
    /* Write all ghost metadata */
    if (!t8_cmesh_load_ghosts (cmesh, fp)) {
      t8_errorf ("Error when opening file %s.\n", filename);
      t8_cmesh_destroy (&cmesh);
      return NULL;
    }
  }
  t8_cmesh_trees_finish_part (cmesh->trees, 0);
  if (!t8_cmesh_load_tree_attributes (cmesh, fp)) {
    t8_errorf ("Error when opening file %s.\n", filename);
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }
  if (cmesh->set_partition) {
    /* Write all ghost metadata */
    if (!t8_cmesh_load_ghost_attributes (cmesh, fp)) {
      t8_errorf ("Error when opening file %s.\n", filename);
      t8_cmesh_destroy (&cmesh);
      return NULL;
    }
  }
  /* Close the file */
  fclose (fp);

  cmesh->committed = 1;
  mpiret = sc_MPI_Comm_rank (comm, &cmesh->mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &cmesh->mpisize);
  SC_CHECK_MPI (mpiret);
  t8_stash_destroy (&cmesh->stash);
  return cmesh;
}

/* Query whether a given process will open a cmesh saved file.
 * This depends on the number of processes, the number of files and
 * the load mode.
 * In load mode simple, if rank i is smaller than the number of files,
 * rank i opens file number i.
 * In load mode bgq, there must be at leas num_files compute nodes
 * and on the i-th node exactly one process opens file i.
 * This function returns true, if the given rank will load a file and
 * it stores the number of the file in file_to_load.
 */
static int
t8_cmesh_load_proc_loads (const int mpirank, const int mpisize, const int num_files, sc_MPI_Comm comm,
                          const t8_load_mode_t mode, int *file_to_load, int num_procs_per_node)
{
  sc_MPI_Comm inter = sc_MPI_COMM_NULL, intra = sc_MPI_COMM_NULL;
  int mpiret, interrank, intrarank, intersize;

  if (num_procs_per_node <= 0 && mode == T8_LOAD_STRIDE) {
    t8_global_infof ("number of processes per node set to 16\n");
    num_procs_per_node = 16;
  }
  /* Fill with invalid value */
  *file_to_load = -1;
  switch (mode) {
  case T8_LOAD_SIMPLE:
    /* In simple mode, process i opens file i, if i < num_files */
    if (mpirank < num_files) {
      /* Open file mpirank */
      *file_to_load = mpirank;
      return 1;
    }
    else {
      /* Do not open any files */
      return 0;
    }
    break;
  case T8_LOAD_BGQ:
    /* In bgq mode on each shared memory node is one process that
     * opens a file. */
    /* Compute internode and intranode communicator */
    sc_mpi_comm_attach_node_comms (comm, 0);
    /* Store internode and intranode communicator */
    sc_mpi_comm_get_node_comms (comm, &intra, &inter);
    /* Abort if we could not compute these communicators */
    SC_CHECK_ABORT (intra != sc_MPI_COMM_NULL && inter != sc_MPI_COMM_NULL, "Could not get proper internode "
                                                                            "and intranode communicators.\n");
    mpiret = sc_MPI_Comm_size (inter, &intersize);
    SC_CHECK_MPI (mpiret);
    /* Check if the number of nodes is at least as big as the number of files */
    SC_CHECK_ABORTF (intersize <= num_files,
                     "Must have more compute nodes "
                     "than files. %i nodes and %i fields are given.\n",
                     intersize, num_files);
    mpiret = sc_MPI_Comm_rank (inter, &interrank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (intra, &intrarank);
    SC_CHECK_MPI (mpiret);
    /* If the current node number i is smaller than the number of files
     * than the first process on this node loads the file i */
    if (interrank < num_files && intrarank == 0) {
      *file_to_load = interrank;
      return 1;
    }
    else {
      return 0;
    }
    break;
  case T8_LOAD_STRIDE:
    /* In Juqueen mode, every 16-th process loads a file. The user should
     * control, that these processes reside on different compute nodes to
     * gain maximal efficiency. */
    SC_CHECK_ABORT (ceil (mpisize / (double) num_procs_per_node) >= num_files,
                    "Too many files for too few processes.\n");
    if (mpirank % num_procs_per_node == 0 && mpirank / num_procs_per_node < num_files) {
      *file_to_load = mpirank / num_procs_per_node;
      return 1;
    }
    else {
      return 0;
    }
  default:
    SC_ABORT_NOT_REACHED ();
  }
  SC_ABORT_NOT_REACHED ();
}

/* After we loaded a cmesh on a part of the processes, we have to correctly
 * set the first tree on the other, empty, cmeshes.
 * In load mode simple, the first tree is just the global number of trees.
 * In Load mode BGQ and JUQUEEN, we have to figure out the first tree of
 * the next bigger nonloading process.
 */
static int
t8_cmesh_load_bigger_nonloading (const int mpirank, const int mpisize, const int num_files, const t8_load_mode_t mode,
                                 sc_MPI_Comm comm, const int num_procs_per_node)
{
  int next_bigger_nonloading;
  sc_MPI_Comm inter = sc_MPI_COMM_NULL, intra = sc_MPI_COMM_NULL;
  sc_MPI_Group intragroup, commgroup;
  int mpiret, interrank, intrarank, intrasize, commrank;
  int rankzero;

  switch (mode) {
  case T8_LOAD_SIMPLE:
    /* In simple mode, the first num_files processes load the cmesh and
     * the rest is empty, this the next bigger nonloading rank is always
     * rank mpisize. */
    next_bigger_nonloading = mpisize;
    break;
  case T8_LOAD_BGQ:
    /* In BGQ mode, on the compute nodes with node-id bigger num_files,
     * the first rank on a node opens the file. */
    /* Get the inter and intra comms and our rank within these comms. */
    sc_mpi_comm_get_node_comms (comm, &intra, &inter);
    mpiret = sc_MPI_Comm_rank (inter, &interrank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (intra, &intrarank);
    SC_CHECK_MPI (mpiret);

    if (interrank >= num_files - 1) {
      /* We are beyond all processes that loaded */
      next_bigger_nonloading = mpisize;
    }
    else {
      /* We need to get the rank in comm of the first process
       * on the next compute node. This is the process with rank 0 in our
       * intra communicator plus the mpisize in the intracommunicator. */
      /* First we get the groups from the communicators */
      mpiret = sc_MPI_Comm_group (intra, &intragroup);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Comm_group (comm, &commgroup);
      SC_CHECK_MPI (mpiret);
      /* We can now transform the rank 0 of the intragroup to get
       * the correct rank in the commgroup */
      rankzero = 0;
      mpiret = sc_MPI_Group_translate_ranks (intragroup, 1, &rankzero, commgroup, &commrank);
      SC_CHECK_MPI (mpiret);
      /* We get the size of the intracommunicator */
      mpiret = sc_MPI_Group_size (intragroup, &intrasize);
      SC_CHECK_MPI (mpiret);
      next_bigger_nonloading = commrank + intrasize;
    }
    break;
  case T8_LOAD_STRIDE:
    /* In Juqueen mode, every 16-th process has opened the file. */
    if (mpirank / num_procs_per_node < num_files - 1) {
      /* If we are in a multiple of 16, where a file was loaded,
       * the first process in the next group did load it. */
      next_bigger_nonloading
        = mpirank - mpirank % num_procs_per_node + num_procs_per_node; /* This is the next number divisible by
                                                                                                   num_procs_per_node */
    }
    else {
      /* Every rank that loaded a file is smaller than ours. */
      next_bigger_nonloading = mpisize;
    }
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return next_bigger_nonloading;
}

/* Load the files fileprefix_0000.cmesh, ... , fileprefix_N.cmesh
 * (N = num_files - 1)
 * on N processes and repartition the cmesh to all calling processes.
 * If N = 1, the cmesh is broadcasted and not partitioned.
 */
t8_cmesh_t
t8_cmesh_load_and_distribute (const char *fileprefix, const int num_files, sc_MPI_Comm comm, const t8_load_mode_t mode,
                              const int procs_per_node)
{
  t8_cmesh_t cmesh;
  char buffer[BUFSIZ];
  int mpiret, mpirank, mpisize;
  int file_to_load;
  int next_bigger_nonloading;
  int did_load;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (mpisize >= num_files);

  /* Try to set the comm type */
  SC_CHECK_ABORT (t8_shmem_init (comm) > 0, "Error in shared memory setup. Could not load cmesh.");
  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);

  /* Use cmesh_bcast, if only one process loads the cmesh: */
  if (num_files == 1) {
    cmesh = NULL;
    if (mpirank == 0) {
      snprintf (buffer, BUFSIZ, "%s_%04d.cmesh", fileprefix, 0);
      cmesh = t8_cmesh_load (buffer, comm);
    }
    cmesh = t8_cmesh_bcast (cmesh, 0, comm);
  }
  else {
    /* More than one process loads a file */

    if (t8_cmesh_load_proc_loads (mpirank, mpisize, num_files, comm, mode, &file_to_load, procs_per_node)) {
      T8_ASSERT (fileprefix != NULL);
      T8_ASSERT (0 <= file_to_load && file_to_load < num_files);
      snprintf (buffer, BUFSIZ, "%s_%04d.cmesh", fileprefix, file_to_load);
      t8_infof ("Opening file %s\n", buffer);
      cmesh = t8_cmesh_load (buffer, comm);
      if (num_files == mpisize) {
        /* Each process has loaded the cmesh and we can return */
        return cmesh;
      }
      did_load = 1;
    }
    else {
      did_load = 0;
      /* On the processes that do not load the cmesh, initialize it
       * with zero local trees and ghosts */
      t8_cmesh_init (&cmesh);
      t8_cmesh_trees_init (&cmesh->trees, 0, 0, 0);
      cmesh->mpirank = mpirank;
      cmesh->mpisize = mpisize;
      t8_stash_destroy (&cmesh->stash);
      /* There are no faces, so we know all about them */
      cmesh->face_knowledge = 3;
      /* There are no tree, thus the first tree is not shared */
      cmesh->first_tree_shared = 0;
      cmesh->committed = 1;
      cmesh->set_partition = 1;
      cmesh->num_local_trees = 0;
    }
    /* The cmeshes on the processes that did not load have to
     * know the global number of trees */
    sc_MPI_Bcast (&cmesh->num_trees, 1, T8_MPI_GLOIDX, 0, comm);
    /* And the dimension */
    sc_MPI_Bcast (&cmesh->dimension, 1, sc_MPI_INT, 0, comm);
    t8_cmesh_gather_trees_per_eclass (cmesh, comm);
    T8_ASSERT (t8_cmesh_is_committed (cmesh));
    /* We now create the cmeshs offset in order to properly
     * set the first tree for the empty processes */
    t8_cmesh_gather_treecount (cmesh, comm);
    if (!did_load) {
      /* Calculate the next bigger nonloading rank. */
      next_bigger_nonloading
        = t8_cmesh_load_bigger_nonloading (mpirank, mpisize, num_files, mode, comm, procs_per_node);
      /* Set the first tree of this process to the first tree of the next nonloading
       * rank */
      cmesh->first_tree
        = t8_offset_first (next_bigger_nonloading, t8_shmem_array_get_gloidx_array (cmesh->tree_offsets));
    }
    /* Since we changed the first tree on some processes, we have to
     * regather the first trees on each process */
    t8_cmesh_gather_treecount (cmesh, comm);
  }
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  return cmesh;
}
