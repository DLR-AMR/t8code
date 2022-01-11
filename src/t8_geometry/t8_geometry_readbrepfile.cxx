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

/* In this file we read occ brep files to store the geometries in an occ geometry. */

#include <t8_eclass.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>
#include <t8_geometry/t8_geometry_readbrepfile.h>
#include <t8_cmesh_readmshfile.h>


#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>

template<typename input_iterator, typename lambda_funk>
input_iterator
t8_geometry_find (input_iterator first, input_iterator last, lambda_funk lambda)
{
  for (input_iterator it = first; it != last; ++it)
  {
    if (lambda(*it))
    {
      return it;
    }
  }
  return last;
}

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

t8_geometry_occ_c* 
t8_geometry_from_brep_file (const char *fileprefix, 
                            sc_hash_t *node_table,
                            int dim,
                            double tol,
                            int debugfile)
{
  #if T8_WITH_OCC
  int                               recombined;
  Standard_Real                     first, last;
  t8_msh_file_node_parametric_t    *cur_node;
  sc_list_t                        *list;
  std::string                       current_file(fileprefix);       
  t8_geometry_occ                  *geom = new t8_geometry_occ (dim, fileprefix, "occ brep derived geometry");
  TopTools_IndexedMapOfShape        occ_shape_vertex_map = geom->t8_geom_get_occ_shape_vertex_map();
  TopTools_IndexedMapOfShape        occ_shape_edge_map = geom->t8_geom_get_occ_shape_edge_map();
  TopTools_IndexedMapOfShape        occ_shape_face_map = geom->t8_geom_get_occ_shape_face_map();
  std::vector<double>               debug_recombined, debug_not_recombined;
  std::vector<std::pair<int, int>>  points_processed, curves_processed, surfaces_processed;

  /* To find the corresponding occ geometry for each parametric msh node 
   * we check each node and insert its parameters into each occ geometry of the right dimension.
   * If the coordinates of the node get returned we know it is the right occ geometry. We save the already found
   * occ geometries so we don't have to process them twice. We also save the indices of each recombined
   * node to change the entity_tag to the position of the corresponding occ geometry in the cmesh geometry. */
  
  /* Iterate over each cmesh node */
  for (size_t hval = 0; hval < node_table->slots->elem_count; ++hval)
  {
    list = (sc_list_t *) sc_array_index (node_table->slots, hval);
    for (sc_link_t *item = list->first; item != NULL; item = item->next)
    {
      cur_node = (t8_msh_file_node_parametric_t *)item->data;
      /* Check if node is parametric or on dim 0. Dim 0 nodes have no parameters and gmsh therefore labels them as not parametric. */
      if (cur_node->parametric || cur_node->entity_dim == 0)
      {
        /* Skip node if the corresponding occ point was already found */
        if (cur_node->entity_dim == 0)
        {
          /* Search for entity_tag in the already processed points */
          std::vector<std::pair<int, int>>::iterator 
          it = t8_geometry_find(points_processed.begin(), 
                                points_processed.end(),
                                [&](std::pair<int, int> a){return (a).first == cur_node->entity_tag;});
          if (it != points_processed.end())
          {
            /* If the entity_tag was already processed and the index of the corresponding point is not 0 we change the 
             * entity_tag of the node to the index of the corresponding point. If the index is 0 we set the node to not parametric. */
            if (it->second > 0)
            {
              cur_node->entity_tag = it->second;
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
              cur_node->parametric = 1;
            }
            else
            {
              cur_node->parametric = 0;
              debug_not_recombined.insert(debug_not_recombined.end(),
                                          cur_node->coordinates,
                                          cur_node->coordinates + 3);
              t8_global_errorf("WARNING: Failed to recombine node with index %i\n", cur_node->index);
            }
            /* Skip node */
            continue;
          }
        }
        /* Skip node if the corresponding occ curve has already been found */
        if (cur_node->entity_dim == 1)
        {
          /* Search for entity_tag in the already processed curves */
          std::vector<std::pair<int, int>>::iterator 
          it = t8_geometry_find(curves_processed.begin(), 
                                curves_processed.end(),
                                [&](std::pair<int, int> a){return (a).first == cur_node->entity_tag;});
          if (it != curves_processed.end())
          {
            /* If the entity_tag was already processed and the index of the corresponding curve is not 0 we change the 
             * entity_tag of the node to the index of the corresponding curve. If the index is 0 we set the node to not parametric. */
            if (it->second > 0)
            {
              cur_node->entity_tag = it->second;
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
            }
            else
            {
              cur_node->parametric = 0;
              debug_not_recombined.insert(debug_not_recombined.end(),
                                          cur_node->coordinates,
                                          cur_node->coordinates + 3);
              t8_global_errorf("WARNING: Failed to recombine node with index %i\n", cur_node->index);
            }
            /* Skip node */
            continue;
          }
        }
        /* Skip node if the corresponding occ surface has already been found */
        if (cur_node->entity_dim == 2)
        {
          /* Search for entity_tag in the already processed surfaces */
          std::vector<std::pair<int, int>>::iterator 
          it = t8_geometry_find(surfaces_processed.begin(), 
                                surfaces_processed.end(),
                                [&](std::pair<int, int> a){return (a).first == cur_node->entity_tag;});
          if (it != surfaces_processed.end())
          {
            /* If the entity_tag was already processed and the index of the corresponding surface is not 0 we change the 
             * entity_tag of the node to the index of the corresponding surface. If the index is 0 we set the node to not parametric. */
            if (it->second > 0)
            {
              cur_node->entity_tag = it->second;
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
            }
            else
            {
              cur_node->parametric = 0;
              debug_not_recombined.insert(debug_not_recombined.end(),
                                          cur_node->coordinates,
                                          cur_node->coordinates + 3);
              t8_global_errorf("WARNING: Failed to recombine node with index %i\n", cur_node->index);
            }
            /* Skip node */
            continue;
          }
        }
        
        /* If the entity_tag hasn't already been processed we search for the corresponding occ geometry */
        recombined = 0;
        gp_Pnt cur_pnt;
        if (cur_node->entity_dim == 0)
        {
          /* Iterate over each occ point */
          for (auto point = occ_shape_vertex_map.cbegin(); point != occ_shape_vertex_map.cend(); ++point)
          {
            /* Check if point and node have the same coordinates */
            cur_pnt = BRep_Tool::Pnt(TopoDS::Vertex(*point));
            if (tol > sqrt(pow(cur_pnt.X() - cur_node->coordinates[0], 2)
                  + pow(cur_pnt.Y() - cur_node->coordinates[1], 2)
                  + pow(cur_pnt.Z() - cur_node->coordinates[2], 2)))
            {
              /* Save recombined entity_tag and occ point index. */
              points_processed.push_back(std::make_pair(cur_node->entity_tag, occ_shape_vertex_map.FindIndex(*point)));
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
              cur_node->entity_tag = occ_shape_vertex_map.FindIndex(*point);
              cur_node->parametric = 1;
              recombined = 1;
              break;
            }
          }
          if (!recombined)
          {
            /* Save entity_tag with 0 to mark, that no matching geometry was found.
             * Save node as not parametric. */
            points_processed.push_back(std::make_pair(cur_node->entity_tag, 0));
            cur_node->parametric = 0;
            debug_not_recombined.insert(debug_not_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
            t8_global_errorf("WARNING: Failed to recombine node with index %i \n", cur_node->index);
          }
        }
        if (cur_node->entity_dim == 1)
        {
          /* Iterate over each occ curve */
          for (auto curve = occ_shape_edge_map.cbegin(); curve != occ_shape_edge_map.cend(); ++curve)
          {
            /* Evaluate occ curve with the parameter of the node and check if the node coordinates get returned */
            if(!BRep_Tool::Degenerated(TopoDS::Edge(*curve)))
            {
              BRep_Tool::Curve(TopoDS::Edge(*curve), first, last)->D0(cur_node->parameters[0], cur_pnt);
              if (tol > sqrt(pow(cur_pnt.X() - cur_node->coordinates[0], 2)
                  + pow(cur_pnt.Y() - cur_node->coordinates[1], 2)
                  + pow(cur_pnt.Z() - cur_node->coordinates[2], 2)))
              {
                /* Save recombined entity_tag and occ curve index. */
                curves_processed.push_back(std::make_pair(cur_node->entity_tag, occ_shape_edge_map.FindIndex(*curve)));
                debug_recombined.insert(debug_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
                cur_node->entity_tag = occ_shape_edge_map.FindIndex(*curve);
                recombined = 1;
                break;
              }
            }
          }
          if (!recombined)
          {
            /* Save entity_tag with 0 to mark that no matching geometry was found.
             * Save node as not parametric. */
            curves_processed.push_back(std::make_pair(cur_node->entity_tag, 0));
            cur_node->parametric = 0;
            debug_not_recombined.insert(debug_not_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
            t8_global_errorf("WARNING: Failed to recombine node with index %i \n", cur_node->index);
          }
        }
        if (cur_node->entity_dim == 2)
        {
          /* Iterate over each occ surface */
          for (auto surface = occ_shape_face_map.cbegin(); surface != occ_shape_face_map.cend(); ++surface)
          {
            /* Evaluate occ surface with the parameters of the node and check if the node coordinates get returned. */
            BRep_Tool::Surface(TopoDS::Face(*surface))->D0(cur_node->parameters[0], cur_node->parameters[1], cur_pnt);
            if (tol > sqrt(pow(cur_pnt.X() - cur_node->coordinates[0], 2)
                + pow(cur_pnt.Y() - cur_node->coordinates[1], 2)
                + pow(cur_pnt.Z() - cur_node->coordinates[2], 2)))
            {
              /* Save recombined entity_tag and occ surface index. */
              surfaces_processed.push_back(std::make_pair(cur_node->entity_tag, occ_shape_face_map.FindIndex(*surface)));
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
              cur_node->entity_tag = occ_shape_face_map.FindIndex(*surface);
              recombined = 1;
              break;
            }
          }
          if (!recombined)
          {
            /* Save entity_tag with 0 to mark, that no matching geometry was found.
             * Save node as not parametric. */
            surfaces_processed.push_back(std::make_pair(cur_node->entity_tag, 0));
            cur_node->parametric = 0;
            debug_not_recombined.insert(debug_not_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
            t8_global_errorf("WARNING: Failed to recombine node with index %i \n", cur_node->index);
          }
        }
      }
    }
  }
  
  #ifdef T8_ENABLE_DEBUG
  /* Write out a .geo file which indicates which nodes were successfully recombined. 
   * The file can be opened with gmsh. */
  if (debugfile)
  {
    std::ofstream geofile;
    int green = 0, red;
    geofile.open (current_file + ".geo");
    for (auto it = debug_recombined.begin(); it != debug_recombined.end(); ++it)
    {
      geofile << "Point(" << green << ") = {" << *it << ", " << *(++it) << ", " << *(++it) << "};" << std::endl;
      ++green;
    }
    red = green;
    for (auto it = debug_not_recombined.begin(); it != debug_not_recombined.end(); ++it)
    {
      geofile << "Point(" << red << ") = {" << *it << ", " << *(++it) << ", " << *(++it) << "};" << std::endl;
      ++red;
    }
    if (debug_recombined.size())
    {
      geofile << "Color Green{Point{" << std::endl;
      for (int i_green = 0; i_green != green; ++i_green)
      {
        geofile << i_green; 
        if (i_green != green - 1) geofile << ",";
        geofile  << std::endl;
      }
      geofile << "};}" << std::endl;
    }
    if (debug_not_recombined.size())
    {
      geofile << "Color Red{Point{" << std::endl;
      for (int i_red = green; i_red != red; ++i_red)
      {
        geofile << i_red;
        if (i_red != red - 1) geofile << ",";
        geofile << std::endl;
      }
      geofile << "};}" << std::endl;
    }
    geofile.close();
  }
  #endif /*T8_ENABLE_DEBUG*/

  return (t8_geometry_occ_c *) geom;
  
  #else /* !T8_WITH_OCC */
  SC_ABORTF("OCC not linked");
  #endif /* T8_WITH_OCC */
}

T8_EXTERN_C_END ();
