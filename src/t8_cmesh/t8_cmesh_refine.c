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

/** \file t8_cmesh_refine.c
 *
 * TODO: document this file
 */

/* TODO: could this file be part of cmesh_commit.c? */

#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_partition.h"


t8_locidx_t
t8_cmesh_refine_new_localid (t8_locidx_t parent_id, int child_id,
                             int factor)
{
  return parent_id * factor + child_id;
}

/* TODO: comment */
/* The eclass is the eclass of the parent tree */
void
t8_cmesh_refine_new_neighbors (t8_locidx_t parent_id,
                               t8_eclass_t eclass, int child_id,
                               t8_locidx_t * neighbor_out, int8_t *ttf_out,
                               int factor)
{
  t8_locidx_t         iface, old_neigh;
  int                 num_faces, F, dim, old_face, orient, new_neigh_childid;


  dim = t8_eclass_to_dimension[eclass];
  F = t8_eclass_max_num_faces[dim];
  num_faces = t8_eclass_num_faces[eclass];
  switch (eclass) {
    case T8_ECLASS_TRIANGLE:
      if (child_id == 3) {
        /* The triangle in the middle is connect to child i on face i.
         * The orientation of each edge is 1 */
        for (iface = 0;iface < num_faces;iface++) {
          /* Set the new local id of the face neighbors */
          neighbor_out[iface] = t8_cmesh_refine_new_localid(parent_id, iface,
                                                            factor);
          /* Set the new orientation and face number */
          ttf_out[iface] = 1 * F + iface;
        }
      }
      else {
        int neigh_childidsA[6] = {1,0,0,2,2,1};
        int neigh_childidsB[6] = {2,2,1,1,0,0};
        /* child_id is 0,1,2 */
        T8_ASSERT (0 <= child_id && child_id <= 2);
        /* Along face number child_id we are connected with the middle triangle
         * which has child_id = 3 */
        neighbor_out[child_id] = t8_cmesh_refine_new_localid (parent_id, 3,
                                                              factor);
        ttf_out[child_id] = 1 * F + child_id;
        for (iface = 0;iface < num_faces;iface++) {
          if (iface == child_id) {
            /* We already set this side */
            continue;
          }
          /* Get the old face neighbor, face number and orientation */
          old_neigh = neighbor_out[iface];
          old_face = ttf_out[iface] % F;
          orient = ttf_out[iface] / F;
          /* depending on the child_id and old_face we set the new neighbors
           * child_id */
          if (2 * child_id + old_face <= 2) {
            new_neigh_childid = neigh_childidsA[orient * 3 + old_face];
          }
          else {
            new_neigh_childid = neigh_childidsB[orient * 3 + old_face];
          }
          /* Compute the new neighbors local id from the old one and the new
           * child id */
          neighbor_out[iface] = t8_cmesh_refine_new_localid (old_neigh,
                                                             new_neigh_childid,
                                                             factor);
          /* Set the new orientation and face number are the same as the old
           * ones so we do not need to change ttf_out */
        }
      }
      break;
    default: SC_ABORTF ("Not implemented %s.\n", "yet");
  }
}

/* Compute the coordinates of the i-th child of a tree of a given class
 * from the coordinates of that three.
 *
 * !!! For Triangles and Tets we use Bey's original order !!!
 *
 *          x_2
 *          /\
 *         /2 \
 *        /____\
 *       /\ 3  /\       For Triangles.
 *      /0 \  / 1\
 *     /____\/____\
 *    x_0          x_1
 *
 */
/* TODO: Implement for prisms, Pyramids, points and lines */
void
t8_cmesh_refine_new_coord (const double *coords_in,
                           t8_eclass_t eclass, int child_id,
                           double * coords_out)
{
  int         num_vertices, ivertex, dim, idim;
  int         coord_lookup[8];

  T8_ASSERT (coords_in != NULL);
  T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_QUAD
             || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET);

  num_vertices = t8_eclass_num_vertices[eclass];
  dim = t8_eclass_to_dimension[eclass];
  switch (eclass){
    case T8_ECLASS_QUAD:
    case T8_ECLASS_HEX:
      for (ivertex = 0;ivertex < num_vertices;ivertex++) {
        for (idim = 0;idim < dim;idim++) {
          /* Xout_i = Xin_i,childid   ,where i=ivertex and x_ij = (x_i+x_j)/2 */
          coords_out[dim * ivertex + idim] =
              (coords_in[dim * ivertex + idim] +
              coords_in[dim * child_id + idim])/2;
        }
      }
      break;
    case T8_ECLASS_TET:
    case T8_ECLASS_TRIANGLE:
      switch (child_id){
        case 0:
        case 1:
        case 2:
        case 3: /* Xout_i = Xin_i,childid   ,where i=ivertex and x_ij = (x_i+x_j)/2 */
          coord_lookup[0] = 0;
          coord_lookup[1] = child_id;
          coord_lookup[2] = 1;
          coord_lookup[3] = child_id;
          coord_lookup[4] = 2;
          coord_lookup[5] = child_id;
          coord_lookup[6] = 3;
          coord_lookup[7] = child_id;
          if (eclass == T8_ECLASS_TRIANGLE && child_id == 3) {
            coord_lookup[0] = 0;
            coord_lookup[1] = 1;
            coord_lookup[2] = 0;
            coord_lookup[3] = 2;
            coord_lookup[4] = 1;
            coord_lookup[5] = 2;
          }
          break;
        case 4: coord_lookup[0] = 0;
          coord_lookup[1] = 1;
          coord_lookup[2] = 0;
          coord_lookup[3] = 2;
          coord_lookup[4] = 0;
          coord_lookup[5] = 3;
          coord_lookup[6] = 1;
          coord_lookup[7] = 3;
          break;
        case 5: coord_lookup[0] = 0;
          coord_lookup[1] = 1;
          coord_lookup[2] = 0;
          coord_lookup[3] = 2;
          coord_lookup[4] = 1;
          coord_lookup[5] = 2;
          coord_lookup[6] = 1;
          coord_lookup[7] = 3;
          break;
        case 6: coord_lookup[0] = 0;
          coord_lookup[1] = 2;
          coord_lookup[2] = 0;
          coord_lookup[3] = 3;
          coord_lookup[4] = 1;
          coord_lookup[5] = 3;
          coord_lookup[6] = 2;
          coord_lookup[7] = 3;
          break;
        case 7: coord_lookup[0] = 0;
          coord_lookup[1] = 2;
          coord_lookup[2] = 1;
          coord_lookup[3] = 2;
          coord_lookup[4] = 1;
          coord_lookup[5] = 3;
          coord_lookup[6] = 2;
          coord_lookup[7] = 3;
          break;
      }

      for (ivertex = 0;ivertex < num_vertices;ivertex++) {
        for (idim = 0;idim < 3;idim++) {
          coords_out[3 * ivertex + idim] =
              (coords_in[3 * coord_lookup[2*ivertex] + idim] +
              coords_in[3 * coord_lookup[2*ivertex + 1] + idim])/2;
        }
      }
      break;
    default: SC_ABORT_NOT_REACHED ();
  }
}

/* Set the number of attributes and the class for each child of a given tree in cmesh_from */
void
t8_cmesh_refine_inittree (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                          t8_locidx_t treeid, t8_locidx_t firstnewtree,
                          int factor)
{
  t8_ctree_t          tree, newtree;
  t8_locidx_t         itree;
  t8_locidx_t        *tree_neighbors;
  int8_t             *ttf;
  size_t              newattr_size;


  tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, treeid,
                                      &tree_neighbors, &ttf);
  newattr_size = t8_cmesh_trees_attribute_size (tree);
  for (itree = 0;itree < factor;itree++) {
    newtree = t8_cmesh_trees_get_tree (cmesh->trees, firstnewtree + itree);
    newtree->eclass = tree->eclass; /* TODO: For pyramid support we will
                                        need to change the eclass here to tets for some trees*/
    newtree->num_attributes = tree->num_attributes;
    newtree->treeid = firstnewtree + itree;
    t8_cmesh_trees_init_attributes (cmesh->trees, firstnewtree + itree,
                                    newtree->num_attributes, newattr_size);
  }
}

/* Set the attributes and face_neighbors for each child of a given tree (treeid) in cmesh_from */
/* The attributes are just copies of the parent's attributes,
 * except for the coordinate attribute if it is set (t8_package_id and key=0).
 * The new coordinates are then calculated.
 */
void
t8_cmesh_refine_tree (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                      t8_locidx_t treeid, t8_locidx_t firstnewtree, int factor)
{
  t8_ctree_t          tree, newtree;
  t8_locidx_t         itree;
  t8_locidx_t        *tree_neighbors, *ntree_neighbors;
  int8_t             *ttf, *nttf;
  size_t              iatt;
  double             *coords;
  t8_attribute_info_struct_t *attr_info;
  t8_stash_attribute_struct_t attr_struct;

  tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, treeid,
                                      &tree_neighbors, &ttf);
  for (itree = 0;itree < factor;itree++) {
    newtree = t8_cmesh_trees_get_tree_ext (cmesh->trees, firstnewtree + itree,
                                           &ntree_neighbors, &nttf);
    /* Set all attributes of the child tree */
    for (iatt = 0;iatt < tree->num_attributes;iatt++) {
      attr_info = T8_TREE_ATTR_INFO (tree, iatt);
      attr_struct.attr_data = T8_TREE_ATTR (tree, attr_info);
      attr_struct.attr_size = attr_info->attribute_size;
      attr_struct.id = cmesh->first_tree + firstnewtree + itree;
      attr_struct.key = attr_info->key;
      attr_struct.package_id = attr_info->package_id;
      attr_struct.is_owned = 0; /* Make sure that the attribute data is copied
                                  from the old tree */
      t8_cmesh_trees_add_attribute (cmesh->trees, 0, &attr_struct,
                                    firstnewtree + itree, iatt);
      if (attr_struct.package_id == t8_get_package_id ()
          && attr_struct.key == 0) {
        /* The attribute was the tree's coordinates.
         * In this case we compute the new coordinates */
        coords = (double *) T8_TREE_ATTR (newtree, T8_TREE_ATTR_INFO (newtree,
                                                                      iatt));
        /* coords is the attribute of newtree */
        t8_cmesh_refine_new_coord (attr_struct.attr_data, tree->eclass, itree,
                                   coords);
        attr_struct.is_owned = 1;
      }
    }
    /* Set all face_neighbors of the child tree */
    t8_cmesh_refine_new_neighbors (treeid, tree->eclass, itree, ntree_neighbors,
                                   nttf, factor);
  }
}

/* Set the class for each child of a given ghost in cmesh_from */
void
t8_cmesh_refine_initghost (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                           t8_locidx_t ghostid, t8_locidx_t firstnewghost,
                           int factor)
{
  t8_cghost_t         ghost, newghost;
  t8_locidx_t         ighost;
  t8_gloidx_t        *ghost_neighbors;
  int8_t             *ttf;


  ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees, ghostid,
                                      &ghost_neighbors, &ttf);
  for (ighost = 0;ighost < factor;ighost++) {
    newghost = t8_cmesh_trees_get_ghost (cmesh->trees, firstnewghost + ighost);
    newghost->eclass = ghost->eclass; /* TODO: For pyramid support we will
                                        need to change the eclass here to tets for some ghosts*/
    newghost->treeid = ghost->treeid * factor + ighost;
  }
}


void
t8_cmesh_refine_ghost (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                       t8_locidx_t ghostid, t8_locidx_t firstnewghost,
                       int factor)
{
  t8_cghost_t         newghost;
  t8_gloidx_t        *ghost_neighbors, *nghost_neighbors;
  int8_t             *ttf, *nttf;
  t8_locidx_t         ighost;
  int                 num_faces;

  (void) t8_cmesh_trees_get_ghost_ext (cmesh_from->trees, ghostid,
                                       &ghost_neighbors, &ttf);
  for (ighost = 0;ighost < factor;ighost++) {
    newghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees,
                                             firstnewghost + ighost,
                                             &nghost_neighbors, &nttf);
    /* Set all face_neighbors of the child ghost */
    num_faces = t8_eclass_num_faces[newghost->eclass];
    /* TODO: The copied face neighbors are not the correct face neighbors.
     *       Set the correct ones */
    memcpy (nghost_neighbors, ghost_neighbors,
            num_faces * sizeof (t8_gloidx_t));
    memcpy (nttf, ttf, num_faces * sizeof (int8_t));

  }
}

void            t8_cmesh_refine (t8_cmesh_t cmesh)
{
  t8_cmesh_t          cmesh_from;
  t8_locidx_t         itree, firstnewtree;
  t8_locidx_t         ighost, firstnewghost;
  int                 dim, factor, factor_ghosts, level, iclass;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->from_method == T8_CMESH_FROM_REFINE);
  T8_ASSERT (cmesh->set_from->committed);
  T8_ASSERT (cmesh->set_from->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0);
  T8_ASSERT (cmesh->set_level == 1); /* levels bigger than 1 are not yet implemented */

  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  dim = cmesh_from->dimension;
  level = cmesh->set_level;
  /* The number of new trees per old tree
   * dim     factor (level = 1)
   *  0         1   (points)
   *  1         2   (lines)
   *  2         4   (quads and triangles)
   *  3         8   (Hexes, prisms, Tets)
   */
  factor = 1 << (dim * level);
  cmesh->num_local_trees = cmesh_from->num_local_trees * factor;
  cmesh->num_trees = cmesh_from->num_trees * factor;
  for (iclass = T8_ECLASS_FIRST;iclass < T8_ECLASS_LAST;iclass++) {
    cmesh->num_trees_per_eclass[iclass] =
        cmesh_from->num_trees_per_eclass[iclass] * factor;
  }
  /* Since we only consider face-ghosts, the numer of new ghosts per old ghosts is the number of
   * face-children of a face. */
  /* The number of new ghosts per old ghosts
   * dim     factor (level = 1)
   *  0         0   (points -> No ghosts)
   *  1         1   (lines -> boundaries are points)
   *  2         2   (quads and triangles -> boundaries are lines)
   *  3         4   (Hexes,prisms and Tets -> boundaries are quads/triangles)
   */

  factor_ghosts = 1 << ((dim - 1) * level); /* The number of new ghosts per old ghosts */
  cmesh->num_ghosts = cmesh_from->num_ghosts * factor_ghosts;
  cmesh->first_tree = cmesh_from->first_tree * factor;
  /* Check for locidx overflow */
  T8_ASSERT ((t8_gloidx_t) cmesh_from->num_local_trees * factor ==
             cmesh->num_local_trees);
  /************************/
  /* Create the new trees */
  /************************/
  /* We create the new cmesh with only one part, independent of the number of
   * parts of cmesh_from */
  t8_cmesh_trees_init (&cmesh->trees, 1, cmesh->num_local_trees,
                       cmesh->num_ghosts);
  t8_cmesh_trees_start_part (cmesh->trees, 0, 0, cmesh->num_local_trees, 0,
                             cmesh->num_ghosts);
  /* Loop over all trees in cmesh_from and set the classes for their children */
  for (itree = 0, firstnewtree = 0;itree < cmesh_from->num_local_trees;
       itree++, firstnewtree += factor) {
    t8_cmesh_refine_inittree (cmesh, cmesh_from, itree, firstnewtree, factor);
  }
  /* Allocate face neihbors and attributes for new trees */
  t8_cmesh_trees_finish_part (cmesh->trees, 0);
  /* Loop over all trees in cmesh_from and set the new face-neighbors and attributes */
  for (itree = 0, firstnewtree = 0;itree < cmesh_from->num_local_trees;
       itree++, firstnewtree += factor) {
    t8_cmesh_refine_tree (cmesh, cmesh_from, itree, firstnewtree, factor);
  }
  for (ighost = 0, firstnewghost = 0;ighost < cmesh_from->num_ghosts;ighost++,
       firstnewghost += factor_ghosts) {
    t8_cmesh_refine_ghost (cmesh, cmesh_from, ighost, firstnewghost,
                           factor_ghosts);
  }
}
