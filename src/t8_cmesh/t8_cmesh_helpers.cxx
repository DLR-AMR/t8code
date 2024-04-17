/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

/** \file t8_cmesh_helpers.c
 *
 * Collection of cmesh helper routines.
 */

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <bits/stdc++.h>
#include <map>

T8_EXTERN_C_BEGIN ();

void
t8_cmesh_set_join_by_vertices (t8_cmesh_t cmesh, const int ntrees, const t8_eclass_t *eclasses, const double *vertices,
                               int **connectivity, const int do_both_directions)
{
  /* If `connectivity` is NULL then the following array gets freed at the end of this routine. */
  int *conn = T8_ALLOC (int, ntrees *T8_ECLASS_MAX_FACES * 3);
  for (int i = 0; i < ntrees * T8_ECLASS_MAX_FACES * 3; i++) {
    conn[i] = -1;
  }

  double min_coord = vertices[0];
  double max_coord = vertices[0];

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];
    const int nverts = t8_eclass_num_vertices[eclass];
    const int edim = t8_eclass_to_dimension[eclass];

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < edim; icoord++) {
        const double coord = vertices[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)];

        if (coord < min_coord) {
          min_coord = coord;
        }

        if (coord > max_coord) {
          max_coord = coord;
        }
      }
    }
  }

  std::multimap<unsigned int, std::pair<int,int>> faces;

  const double eps = 1e-9;

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];
    const int nverts = t8_eclass_num_vertices[eclass];

    /* Get the number of faces of this element. */
    const int nfaces = t8_eclass_num_faces[eclass];

    /* Loop over all faces of the current cmesh element. */
    for (int iface = 0; iface < nfaces; iface++) {
      /* Get the number of vertices per face of this element. */
      const int nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[eclass][iface]];

      unsigned int uintface = 0;

      /* Loop over the vertices of the current element's face. */
      for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
        /* Map from a face vertex id to the element vertex id. */
        const int ivert = t8_face_vertex_to_tree_vertex[eclass][iface][iface_vert];

        for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
          const int index = T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord);
          const double rescaled = (vertices[index] - min_coord) / (max_coord - min_coord) / eps;
          const unsigned int uintvert = static_cast<unsigned int>(rescaled);

          /* Simple hash function. */
          uintface = uintface ^ (uintvert << icoord);
        }
      }

      // t8_global_productionf ("itree = %d | iface = %d | uintface = %u | count = %ld\n", itree, iface, uintface, faces.count(uintface));

      auto range = faces.equal_range(uintface);

      for (auto it = range.first; it != range.second; ++it) {

        const int neigh_itree = std::get<0>(it->second);
        const int neigh_iface = std::get<1>(it->second);

        // t8_global_productionf ("itree = %d | iface = %d | uintface = %u\n", neigh_itree, neigh_iface, it->first);

        /* If we already checked the two faces, we can
         * skip the computations here since we only need the connectivity in one direction. */
        if (!do_both_directions
            && conn[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_FACES, 3, neigh_itree, neigh_iface, 0)] > -1) {
          continue;
        }

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
              if (fabs (face_vert - neigh_face_vert) < 10.0 * T8_PRECISION_EPS) {
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

          // t8_global_productionf ("itree = (%d, %d) | iface = (%d, %d) | uintface = %u\n", itree, neigh_itree, iface, neigh_iface, uintface);

          /* Compute the orientation of the face-to-face connection.
           * Face corner 0 of the face with the lower face direction connects
           * to a corner of the other face. The number of this corner is the
           * orientation code. */
          int orientation = -1;
          int smaller_bigger_face_condition = -1;

          int compare = t8_eclass_compare (eclass, neigh_eclass);
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

          break;
        }
      } /* Loop over faces with identical hash. */

      faces.insert(std::make_pair(uintface, std::make_pair(itree,iface)));

      // t8_global_productionf ("itree = (%d, %d) | iface = (%d, %d) | uintface = %u\n", faces.itree, iface, uintface);
      // t8_global_productionf ("itree = %d | iface = %d | uintface = %u\n", itree, iface, uintface);
    } /* Loop over faces. */
  } /* Loop over trees. */

  /* Transfer the computed face connectivity to the `cmesh` object. */
  if (cmesh != NULL) {
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
  if (connectivity == NULL) {
    T8_FREE (conn);
  }
  else {
    *connectivity = conn;
  }
}

T8_EXTERN_C_END ();
