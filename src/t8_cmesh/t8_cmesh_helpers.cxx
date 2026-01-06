/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023, 2024 the developers

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

/** \file t8_cmesh_helpers.cxx
 *
 * Collection of cmesh helper routines.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <vector>
#include <map>

void
t8_cmesh_set_join_by_vertices (t8_cmesh_t cmesh, const t8_gloidx_t ntrees, const t8_eclass_t *eclasses,
                               const double *vertices, int **connectivity, const int do_both_directions)
{
  /* If `connectivity` is NULL then the following array gets freed at the end of this routine. */
  int *conn = T8_ALLOC (int, ntrees *T8_ECLASS_MAX_FACES * 3);
  for (int i = 0; i < ntrees * T8_ECLASS_MAX_FACES * 3; i++) {
    conn[i] = -1;
  }

  /* Compute minimum and maximum of the cmesh domain. */
  double min_coord = vertices[0];
  double max_coord = vertices[0];

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];
    const int nverts = t8_eclass_num_vertices[eclass];
    const int edim = t8_eclass_to_dimension[eclass];

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < edim; icoord++) {
        const double coord
          = vertices[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)];

        if (coord < min_coord) {
          min_coord = coord;
        }

        if (coord > max_coord) {
          max_coord = coord;
        }
      }
    }
  }

  /* Scaled tolerance to decide if two doubles are equal. */
  const double tolerance = 10.0 * T8_PRECISION_EPS * std::abs (max_coord - min_coord);

  /* Setup hash table `faces` mapping a hash key to a pair containing `(itree, iface)`. */
  std::multimap<unsigned long, std::pair<int, int>> faces;

  /* `num_bins` should be more than enough for (almost) all cases.
   * I.e., 2^P4EST_QMAXLEVEL =~ 1.073e9.
   */
  const double num_bins = 1e9; /* Number of bins. */
  const double inverse_bin_size = num_bins / (max_coord - min_coord);

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];

    /* Get the number of faces of this element. */
    const int nfaces = t8_eclass_num_faces[eclass];

    /* Loop over all faces of the current cmesh element. */
    for (int iface = 0; iface < nfaces; iface++) {
      /* Get the number of vertices per face of this element. */
      const int nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[eclass][iface]];

      /* Compute the hash key. The idea is to convert the rescaled vertices of a tree face to
       * long integers and add them up. We apply a bit of seeding by also adding `icoord`. See below. */
      unsigned long hash = 0;

      /* Loop over the vertices of the current element's face. */
      for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
        /* Map from a face vertex id to the element vertex id. */
        const int ivert = t8_face_vertex_to_tree_vertex[eclass][iface][iface_vert];

        for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
          const int index = T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord);
          const double rescaled = (vertices[index] - min_coord) * inverse_bin_size;

          /* Simple hash function. */
          hash = hash + icoord + static_cast<unsigned long> (rescaled + 0.5);
        }
      }

      /* Loop over all pre-registered faces with the same hash. */
      auto range = faces.equal_range (hash);
      for (auto it = range.first; it != range.second; ++it) {
        /* Query potentially neighboring `itree` and `iface`. */
        const int neigh_itree = std::get<0> (it->second);
        const int neigh_iface = std::get<1> (it->second);

        /* Retrieve the potentially neighboring element class. */
        const t8_eclass_t neigh_eclass = eclasses[neigh_itree];

        /* Get the number of vertices per face of potentially neighboring element. */
        const int neigh_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[neigh_eclass][neigh_iface]];

        /* If the number of face vertices do not match we can skip. */
        if (nface_verts != neigh_nface_verts) {
          continue;
        }

        /* The order of the encountered face vertices is needed for computing
         * the orientation later on. Prepare the array for that here. */
        int face_vert_order[T8_ECLASS_MAX_EDGES_2D];
        for (int i = 0; i < T8_ECLASS_MAX_EDGES_2D; i++) {
          face_vert_order[i] = -1;
        }

        int match_count = 0; /* This tracks the number of matching vertices. */
        /* Loop over the vertices of the current element's face. */
        for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
          /* Map from a face vertex id to the element vertex id. */
          const int ivert = t8_face_vertex_to_tree_vertex[eclass][iface][iface_vert];

          /* Loop over the vertices of the potentially neighboring element's face. */
          for (int neigh_iface_vert = 0; neigh_iface_vert < neigh_nface_verts; neigh_iface_vert++) {
            /* Map from a face vertex id to the element vertex id. */
            const int neigh_ivert = t8_face_vertex_to_tree_vertex[neigh_eclass][neigh_iface][neigh_iface_vert];

            int match_count_per_coord = 0; /* Tracks the matching of x, y and z coordinates of two vertices. */
            for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
              /* Retrieve the x, y or z component of the face vertex
               * coordinate from the vertices array. */

              const double face_vert
                = vertices[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)];
              const double neigh_face_vert = vertices[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM,
                                                                   neigh_itree, neigh_ivert, icoord)];

              /* Compare the coordinates with some tolerance. */
              if (std::abs (face_vert - neigh_face_vert) < tolerance) {
                match_count_per_coord++;
              }
            }

            /* In case all x, y and z components match we increase the match_count variable. */
            if (match_count_per_coord == T8_ECLASS_MAX_DIM) {
              match_count++;
              /* Store the encountered face vertex order for later use. */
              face_vert_order[iface_vert] = neigh_iface_vert;
              continue;
            }
          }
        }

        /* If the number of matching face vertices is equal to the actual number of the face's vertices
         * we interpret this as a face-to-face connection between two elements. */
        if (match_count == nface_verts) {
          /* Compute the orientation of the face-to-face connection.
           * Face corner 0 of the face with the lower face direction connects
           * to a corner of the other face. The number of this corner is the
           * orientation code. */
          int orientation = -1;
          int smaller_bigger_face_condition = -1;

          const int compare = t8_eclass_compare (eclass, neigh_eclass);
          if (compare < 0) {
            /* This tree class is smaller than neigh. tree class. */
            smaller_bigger_face_condition = 1;
          }
          else if (compare > 0) {
            /* This tree class is bigger than neigh. tree class. */
            smaller_bigger_face_condition = 0;
          }
          else {
            /* This tree class is the same as the neigh. tree class. 
               Then the face with the smaller face id is the smaller one. */
            smaller_bigger_face_condition = iface < neigh_iface;
          }

          if (smaller_bigger_face_condition) {
            orientation = face_vert_order[0];
          }
          else {
            for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
              if (0 == face_vert_order[iface_vert]) {
                orientation = iface_vert;
                break;
              }
            }
          }

          /* Store the results. */
          conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 0)] = neigh_itree;
          conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 1)] = neigh_iface;
          conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 2)] = orientation;

          if (do_both_directions) {
            conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, neigh_itree, neigh_iface, 0)] = itree;
            conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, neigh_itree, neigh_iface, 1)] = iface;
            conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, neigh_itree, neigh_iface, 2)] = orientation;
          }

          break;
        }
      } /* Loop over faces with identical hash. */

      /* Register the current pair of `itree` and `iface` with given `hash` in the hash table. */
      faces.insert (std::make_pair (hash, std::make_pair (itree, iface)));
    } /* Loop over faces. */
  }   /* Loop over trees. */

  /* Transfer the computed face connectivity to the `cmesh` object. */
  if (cmesh != nullptr) {
    for (int itree = 0; itree < ntrees; itree++) {
      const t8_eclass_t eclass = eclasses[itree];
      const int nfaces = t8_eclass_num_faces[eclass];

      for (int iface = 0; iface < nfaces; iface++) {
        const int neigh_itree = conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 0)];
        const int neigh_iface = conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 1)];
        const int orientation = conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, itree, iface, 2)];

        if (neigh_itree > -1) {
          t8_cmesh_set_join (cmesh, itree, neigh_itree, iface, neigh_iface, orientation);
        }
      }
    }
  }

  /* Pass the `conn` array to the caller if asked for. */
  if (connectivity == nullptr) {
    T8_FREE (conn);
  }
  else {
    *connectivity = conn;
  }
}

void
t8_cmesh_set_join_by_stash (t8_cmesh_t cmesh, int **connectivity, const int do_both_directions)
{
  /* Get some pointers to the cmesh's stash */
  const sc_array_t *classes = &(cmesh->stash->classes);
  const sc_array_t *attributes = &(cmesh->stash->attributes);
  const t8_gloidx_t ntrees = classes->elem_count;
  const size_t num_attributes = attributes->elem_count;
  std::vector<t8_eclass_t> eclasses (ntrees);
  std::vector<double> vertices (ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM, 0.0);

  /* Retrieve the eclasses of the trees from the cmesh's stash */
  for (t8_gloidx_t itree = 0; itree < ntrees; itree++) {
    const t8_stash_class_struct_t *entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx (classes, itree);
    eclasses[entry->id] = entry->eclass;
  }

  /* Retrieve the vertices of the trees from the stash */
  for (size_t iattribute = 0; iattribute < num_attributes; iattribute++) {
    const t8_stash_attribute_struct_t *entry
      = (t8_stash_attribute_struct_t *) t8_sc_array_index_locidx (attributes, iattribute);
    /* If the attribute is a vertex attribute, copy it into a separate vector */
    if (entry->key == T8_CMESH_VERTICES_ATTRIBUTE_KEY && entry->package_id == t8_get_package_id ()) {
      const size_t index = T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, entry->id, 0, 0);
      memcpy (&vertices[index], entry->attr_data, entry->attr_size);
    }
  }

  /* Let t8_cmesh_set_join_by_vertices join the trees */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, eclasses.data (), vertices.data (), connectivity, do_both_directions);
}
