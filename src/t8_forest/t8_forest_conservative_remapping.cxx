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

/** file t8_forest_vtk.h
 */

/* TODO: Document this file */

#ifndef T8_FOREST_CONSERVATIVE_REMAPPING
#define T8_FOREST_CONSERVATIVE_REMAPPING

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h> 
#include <vector>

T8_EXTERN_C_BEGIN ();

/* Search query, a set of points. */
typedef struct
{ 
  std::vector<double[3]> coordinates;/* The corners of our cell. */
  int num_corners;
  std::vector<t8_locidx_t> *intersection_cell_indices;
  std::vector<t8_locidx_t> *intersection_cell_treeid;
  double volume_cell;
  double volume_intersection_cells;

} t8_cell_corners_t;

/* Additional user data that we process during search.
 * For each element we count the number of particles that it contains
 * and we count the total number of elements that we constructed during search. */
typedef struct
{
  t8_cell_corners_t *cornerns_per_cell;   /* List of each cell including the corners. */
  //Liste (Anzahl Zellen Zielgitter) von dynamischen structuren fuer die Zell-IDs des Ursprungsgitters 
  
} t8_search_user_data_t;

/*
 * The search callback.
 * It will be called once per element and generally decides whether or not 
 * to continue the search with the children of the element.
 * Since we will continue as long as there are corners left to consider,
 * we always return 1 here.
 * The search will then only stop when no queries are active (thus, no corners
 * could be in this element) or the element is a leaf element.
 */
static int
t8_search_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                             t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index, void *query,
                             size_t query_index)
{
  T8_ASSERT (query == NULL);

  /* Get a pointer to our user data. */
  t8_search_user_data_t *user_data = (t8_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  return 1;
}

/* The query callback. This will be called for each element once per active query
 * (= particle that may be inside this element).
 * The return value determines whether or not the query remains active for the children
 * of this element.
 * In our example this is the case if the particle is inside the element.
 * The additional input parameter 'is_leaf' will be true if the given element is a leaf
 * element (= an actual element in our forest, not a 'virtual' element in the hierarchy).
 * If the element is a leaf and the particle is contained in it, then we will increase
 * a counter for this element by one.
 * These counters are provided in an sc_array as user data of the input forest.
 */
static int
t8_tutorial_search_query_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index,
                                   void *query, size_t query_index)
{
  int corner_is_inside_element;
  t8_cell_corners_t *corners_of_cell = (t8_cell_corners_t *) query;
  //T8_ASSERT (user_data != NULL);
  T8_ASSERT (query != NULL);
  /* Numerical tolerance for the is_inside_element check. */
  const double tolerance = 1e-8;

  for( int icorner = 0; icorner < corners_of_cell->num_corners ; icorner++ )
  {
    corner_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, corners_of_cell->coordinates[icorner], tolerance);
    if (corner_is_inside_element) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to 
        corners_of_cell->intersection_cell_indices->push_back( element_index );
        corners_of_cell->intersection_cell_treeid->push_back( ltreeid );
      }
      /* The particles is inside the element. This query should remain active.
       * If this element is not a leaf the search will continue with its children. */
      return 1;
    }
    /* The particle is not inside the element. Deactivate this query.
     * If no active queries are left, the search will stop for this element and its children. */
    return 0;
  }
}


void //rueckgabewert anpassen oder als pointer uebergeben - Punktwolke
cell_intersection( t8_forest_t forest1, t8_forest_t forest2, t8_element_t *elem1, t8_element_t *elem2 )
{
//t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id, const t8_element_t *element, int corner_number,
//                              double *coordinates)
//ltree_id in t8_cell_corners_t -> uebergeben

//Eckpunkte des Schnittes zurueckgeben

}

/*
 * We are looking at the intersection between two convex cells. 
 * Thus, the result is again a convex shape. We can build a Tet-Mesh out of these
 * points. Therefore, we chose one point and search for the three points, that are closest
 * to this point. If these are not planar, they are the first tet.
 * From each surface we now search for the point, that is closest to this. The triangle and 
 * the closest point again are a tet. We have to test, if the tet intersects another
 * cell. Then it is no tet of the mesh and the tested triangle is an outer surface and has no
 * neighbor.
 */
void 
point_cloud_to_volume( std::vector<double[3]> points )
{
  // triangulation https://www.kiv.zcu.cz/site/documents/verejne/vyzkum/publikace/technicke-zpravy/2002/tr-2002-02.pdf
}

static sc_array *
t8_build_corners( t8_forest_t forest )
{
  int itree, ielem;
  sc_array *corners = sc_array_new_count (sizeof (t8_cell_corners_t), t8_forest_get_global_num_elements( forest ));
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
    }
  }
}

void
t8_forest_conservative_remapping_planar( t8_forest_t forest_old, t8_forest_t forest_new )
{
  int itree, ielem;
  sc_array *corners = t8_build_corners(forest_new);
  //1) Iterieren ueber Zellen des neuen Forests
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      //2) Zellen des alten Forests suchen, in dem die Eckpunkte der jeweiligen Zelle enthalten sind
      t8_forest_search(forest_old, t8_search_callback, t8_tutorial_search_query_callback, corners );
    }
  }
      //3) Schnitt der Zellen bilden

      //4) Volumen vergleichen

        //5) Schnittpunkte von Kanten berechnen

        //6) Zellen des alten Forests suchen, in dem die Schnittpunkte enthalten sind

        //7) Schnitt der Zellen bilden

        //8) Volumen vergleichen

          //9) Eckpunkte der alten Zellen in der neuen Zelle suchen

          //10) Zellen suchen, die diese Eckpunkte als Eckpunkte haben

          //11) Schnitt der Zellen bilden

          //12) Volumen vergleichen -> nicht gleich? -> Fehler
  
    //13) Berechnen des Wertes der neuen Zelle  
}

T8_EXTERN_C_END ();

/* Volumen berechnen */
//

#endif /* !T8_FOREST_CONSERVATIVE_REMAPPING */