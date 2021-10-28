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
#include <iostream>
#include <vector>
#include <array>
#include <optional>
#include <tuple>
#include <algorithm>
#include <utility>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopExp_Explorer.hxx>

TopoDS_Shape
t8_geometry_import_brep (std::string fileprefix)
{
  TopoDS_Shape        shape;
  BRep_Builder        builder;

  std::ifstream is(fileprefix + ".brep");
  BRepTools::Read(shape, is, builder);
  is.close();
  if(shape.IsNull())
  {
    t8_global_errorf("WARNING: Could not read brep file or brep file contains no shape \n");
  }
  return shape;
}

void 
t8_geometry_extract_brep(const TopoDS_Shape &shape, 
                         std::vector<std::optional<Handle_Geom_Curve>> &edges,
                         std::vector<Handle_Geom_Surface> &surfaces)
{
    TopExp_Explorer explorer;
    for (explorer.Init(shape, TopAbs_EDGE); explorer.More(); explorer.Next())
    {
        Standard_Real first, last;
        if(!BRep_Tool::Degenerated(TopoDS::Edge(explorer.Current())))
        {
            edges.push_back({BRep_Tool::Curve(TopoDS::Edge(explorer.Current()), first, last)});
        }
        else
        {
            edges.push_back({});
        }
    }
    for (explorer.Init(shape, TopAbs_FACE); explorer.More(); explorer.Next())
    {
        surfaces.push_back(BRep_Tool::Surface(TopoDS::Face(explorer.Current())));
    }
}

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
                            sc_hash_t * node_table,
                            int dim,
                            int debugfile)
{
  int                                               recombined;
  t8_msh_file_node_parametric_t                    *cur_node;
  sc_list_t                                        *list;               
  t8_geometry_occ                                  *geom = new t8_geometry_occ (dim, "occ brep derived geometry");
  TopoDS_Shape                                      shape;
  std::vector<std::optional<Handle_Geom_Curve>>     curves;
  std::vector<Handle_Geom_Surface>                  surfaces; 
  std::vector<t8_msh_file_node_parametric_t *>      to_edit;
  std::vector<double>                               debug_recombined, debug_not_recombined;
  std::vector<std::pair<std::optional<Handle_Geom_Curve>, int>>    curves_processed;
  std::vector<std::pair<std::optional<Handle_Geom_Surface>, int>>  surfaces_processed;

  std::string current_file(fileprefix);
  shape = t8_geometry_import_brep (current_file);
  t8_geometry_extract_brep(shape, curves, surfaces);
  t8_global_errorf("%s \n", current_file);
  std::cout << current_file << " " << curves.size() << " " << surfaces.size() << std::endl;

  for (size_t hval = 0; hval < node_table->slots->elem_count; ++hval)
  {
    list = (sc_list_t *) sc_array_index (node_table->slots, hval);
    if (list->elem_count > 0)
    {
      cur_node = (t8_msh_file_node_parametric_t *)list->first->data;
      if (cur_node->entity_dim == 1 || cur_node->entity_dim == 2)
      {
        if (cur_node->entity_dim == 1)
        {
          std::vector<std::pair<std::optional<Handle_Geom_Curve>, int>>::iterator 
          it = t8_geometry_find(curves_processed.begin(), 
                                curves_processed.end(),
                                [&](std::pair<std::optional<Handle_Geom_Curve>, int> a){return (a).second == cur_node->entity_tag;});
          if (it != curves_processed.end())
          {
            if (it->first)
            {
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
            }
            else
            {
              debug_not_recombined.insert(debug_not_recombined.end(),
                                          cur_node->coordinates,
                                          cur_node->coordinates + 3);
              t8_global_errorf("WARNING: Failed to recombine node with index %i list\n", cur_node->index);
            }
            continue;
          }
        }
        else
        {
          std::vector<std::pair<std::optional<Handle_Geom_Surface>, int>>::iterator 
          it = t8_geometry_find(surfaces_processed.begin(), 
                                surfaces_processed.end(), 
                                [&](std::pair<std::optional<Handle_Geom_Surface>, int> a){return a.second == cur_node->entity_tag;});
          if (it  != surfaces_processed.end())
          {
            if (it->first)
            {
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
            }
            else
            {
              debug_not_recombined.insert(debug_not_recombined.end(),
                                          cur_node->coordinates,
                                          cur_node->coordinates + 3);
              t8_global_errorf("WARNING: Failed to recombine node with index %i list\n", cur_node->index);
            }
            continue;
          }
        }
        
        recombined = 0;
        gp_Pnt cur_pnt;
        if (cur_node->entity_dim == 1)
        {
          for (std::vector<std::optional<Handle_Geom_Curve>>::iterator it = curves.begin(); it != curves.end(); ++it)
          {
            if (*it)
            {
              it->value()->D0(cur_node->parameters[0], cur_pnt);
              if (1e-3 > abs(cur_pnt.X() - cur_node->coordinates[0]) &&
                  1e-3 > abs(cur_pnt.Y() - cur_node->coordinates[1]) &&
                  1e-3 > abs(cur_pnt.Z() - cur_node->coordinates[2]))
              {
                curves_processed.push_back(std::make_pair(std::optional<Handle_Geom_Curve>{*it}, cur_node->entity_tag));
                curves.erase(it);
                debug_recombined.insert(debug_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
                recombined = 1;
                break;
              }
            }
          }
          if (!recombined)
          {
            curves_processed.push_back(std::make_pair(std::optional<Handle_Geom_Curve>{}, cur_node->entity_tag));
            debug_not_recombined.insert(debug_not_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
          }
        }
        if (cur_node->entity_dim == 2)
        {
          for (std::vector<Handle_Geom_Surface>::iterator it = surfaces.begin(); it != surfaces.end(); ++it)
          {
            (*it)->D0(cur_node->parameters[0], cur_node->parameters[1], cur_pnt);
            if (1e-3 > abs(cur_pnt.X() - cur_node->coordinates[0]) &&
                1e-3 > abs(cur_pnt.Y() - cur_node->coordinates[1]) &&
                1e-3 > abs(cur_pnt.Z() - cur_node->coordinates[2]))
            {
              surfaces_processed.push_back(std::make_pair(std::optional<Handle_Geom_Surface>{*it}, cur_node->entity_tag));
              surfaces.erase(it);
              debug_recombined.insert(debug_recombined.end(),
                                      cur_node->coordinates,
                                      cur_node->coordinates + 3);
              recombined = 1;
              break;
            }
          }
          if (!recombined)
          {
            surfaces_processed.push_back(std::make_pair(std::optional<Handle_Geom_Surface>{}, cur_node->entity_tag));
            debug_not_recombined.insert(debug_not_recombined.end(),
                                        cur_node->coordinates,
                                        cur_node->coordinates + 3);
          }
        }
        if (!recombined)
        {
          t8_global_errorf("WARNING: Failed to recombine node with index %i \n", cur_node->index);
        }
      }
    }
  }

  std::ofstream geofile;
  int green = 0, red;
  geofile.open (current_file + ".geo");
  for (std::vector<double>::iterator it = debug_recombined.begin(); it != debug_recombined.end(); ++it)
  {
    geofile << "Point(" << green << ") = {" << *it << ", " << *(++it) << ", " << *(++it) << "};" << std::endl;
    ++green;
  }
  red = green;
  for (std::vector<double>::iterator it = debug_not_recombined.begin(); it != debug_not_recombined.end(); ++it)
  {
    geofile << "Point(" << red << ") = {" << *it << ", " << *(++it) << ", " << *(++it) << "};" << std::endl;
    ++red;
  }
  geofile << "Color Green{Point{" << std::endl;
  for (int i_green = 0; i_green != green; ++i_green)
  {
    geofile << i_green; 
    if (i_green != green - 1) geofile << ",";
    geofile  << std::endl;
  }
  geofile << "};}" << std::endl;
  geofile << "Color Red{Point{" << std::endl;
  for (int i_red = green; i_red != red; ++i_red)
  {
    geofile << i_red;
    if (i_red != red - 1) geofile << ",";
    geofile << std::endl;
  }
  geofile << "};}" << std::endl;
  geofile.close();

  return (t8_geometry_occ_c *) geom;
}

T8_EXTERN_C_END ();
