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

#include <cmath>
#include <vector>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>

std::array <double, 3> getCoordsUpperRing(int side) {
    if (side < 0 || side > 3) {
        return {};
    }

    const double phi = M_PI / 2;
    std::array <double, 3> coords;

    coords[0] = cos(side * phi);
    coords[1] = sin(side * phi);
    coords[2] = 2.0 / 3;

    return coords;
}

std::array <double, 3> getCoordsBottomRing (int side){
    if (side < 0 || side > 3) {
        return {};
    }
    std::array <double, 3>  coords;
    const double phi = M_PI / 2 ;
    coords[0]= cos(side*phi);
    coords[1]= sin(side*phi);
    coords[2]= -2.0/3;
    return coords;
}


std::array <double, 3>  getCoordsEquator (int side){
    if(side < 0 || side > 3) {
        return {};
    }
    std::array <double, 3>  coords;
    const double phi = M_PI / 2 ;
    const double shift = M_PI / 4;

    coords[0]= cos(side*phi - shift);
    coords[1]= sin(side*phi - shift);
    coords[2]= 0;
    return coords;
}

t8_cmesh_t
t8_cmesh_new_healpix (sc_MPI_Comm comm)
{
    /* Initialization of the mesh */
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);

    // t8_cmesh_register_geometry<t8_geometry_linear> (cmesh); /* Use spherical geometry. */

    const int ntrees = 12;
    const int nverts = 4; /* Number of vertices per cmesh element. */
    std::vector <double> verts;
    t8_eclass_t all_eclasses[ntrees];
    std::vector <double> all_verts;
    t8_cmesh_register_geometry<t8_geometry_healpix> (cmesh);
    /* Defitition of the tree class. */
    for (int itree = 0; itree < ntrees; itree++) {
        t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
        all_eclasses[itree] = T8_ECLASS_QUAD;
    }
    std::array <double, 3>  northPole = {0,0,1};
    std::array <double, 3>  southPole = {0,0,-1};
    int itree = 0;
    // every side section has 4 quads, this builds the upper layer
    // build upper section
for (int  side = 0; side < 4; side++) {
    verts.clear();

    auto equator = getCoordsEquator((side+1) % 4);
    auto upper0  = getCoordsUpperRing(side % 4);
    auto upper1  = getCoordsUpperRing((side + 1) % 4);

    verts.insert(verts.end(), northPole.begin(), northPole.end());
    verts.insert(verts.end(), upper1.begin(), upper1.end());
    verts.insert(verts.end(), upper0.begin(), upper0.end());
    verts.insert(verts.end(), equator.begin(), equator.end());

    all_verts.insert(std::end(all_verts), std::begin(verts), std::end(verts));
    t8_cmesh_set_tree_vertices(cmesh, itree, verts.data(), nverts);

    itree++;
}

// build middle section
for (int side = 0; side < 4; side++) {
    verts.clear();

    auto bottom0  = getCoordsBottomRing(side % 4);
    auto equator1 = getCoordsEquator((side + 1) % 4);
    auto equator0 = getCoordsEquator(side % 4);
    auto upper0   = getCoordsUpperRing(side % 4);

    verts.insert(verts.end(), bottom0.begin(), bottom0.end());
    verts.insert(verts.end(), equator0.begin(), equator0.end());
    verts.insert(verts.end(), equator1.begin(), equator1.end());
    verts.insert(verts.end(), upper0.begin(), upper0.end());

    all_verts.insert(std::end(all_verts), std::begin(verts), std::end(verts));
    t8_cmesh_set_tree_vertices(cmesh, itree, verts.data(), nverts);

    itree++;
}

// build lower section
for (int side = 0; side < 4; side++) {
    verts.clear();

    auto bottom1 = getCoordsBottomRing((side + 1) % 4);
    auto bottom0 = getCoordsBottomRing(side % 4);
    auto equator = getCoordsEquator((side + 1) % 4);

    verts.insert(verts.end(), southPole.begin(), southPole.end());
    verts.insert(verts.end(), bottom0.begin(), bottom0.end());
    verts.insert(verts.end(), bottom1.begin(), bottom1.end());
    verts.insert(verts.end(), equator.begin(), equator.end());

    all_verts.insert(std::end(all_verts), std::begin(verts), std::end(verts));
    t8_cmesh_set_tree_vertices(cmesh, itree, verts.data(), nverts);

    itree++;
}

  t8_cmesh_set_join_by_vertices(cmesh, 12, all_eclasses, all_verts.data(), nullptr, 0);
  T8_FREE(all_eclasses);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
