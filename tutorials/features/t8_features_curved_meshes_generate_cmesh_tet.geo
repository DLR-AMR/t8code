/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* This file will create a brep and msh file of a Joukowsky airfoil.
 * We can then use these files to generate a curved adaptive mesh in
 * t8code.
*/

/* We select the OpenCASCADE geometry kernel because we use the same kernel
 * in t8code. Furthermore, it has a bigger functionality than the inbuilt
 * GEO kernel. */
SetFactory("OpenCASCADE");

/* Now, we can define the points of our airfoil. The coordinates
 * are calculated via the Joukowsky transform. */
lc = 1;
center_chi = -0.15;
center_eta = 0.08;
radius = Sqrt(((1 - center_chi) * (1 - center_chi)) + center_eta);
delta_angle = (2 * Pi) / 100;

Macro XYCoords
  chi = (radius * Cos(angle)) + center_chi;
  eta = (radius * Sin(angle)) + center_eta;
  If (((chi * chi) + (eta * eta)) == 0)
    x = chi;
    y = eta;
  Else
    x = chi * (1 + (1 / ((chi * chi) + (eta * eta))));
    y = eta * (1 - (1 / ((chi * chi) + (eta * eta))));
  EndIf
Return

/* Generate 100 point, defining the airfoil geometry via Joukowsky transform. */
angle =  0;
For i In {0 : 100}
  Call XYCoords;
  Point(1000 + i) = {x, y, 0, lc};
  angle = angle + delta_angle;
EndFor

/* Build a spline fit for the points. Starting and ending with the
 * trailing edge. */
Spline(100) = {1000:1099, 1000};

/* Definition of the corner points of the flow domain. */
Point(2040) = {-3, -2, 0, 1.0};
Point(2041) = {-3, 2, 0, 1.0};
Point(2042) = {4, 2, 0, 1.0};
Point(2043) = {4, -2, 0, 1.0};

/* Build lines through the points. */
Line(101) = {2040, 2041};
Line(102) = {2041, 2042};
Line(103) = {2042, 2043};
Line(104) = {2043, 2040};

/* Close lines to loop */
Curve Loop(1) = {101, 102, 103, 104};
Curve Loop(2) = {100};

/* Make surface with profile cutout */
Plane Surface(1) = {1, 2};
/* And extrude the surface to get a volume. */
Extrude {0, 0, 7} {
  Surface{1};
}
/* Due to the indexing behavior if Gmsh, we have to save the .brep file
 * and reopen it again. */
Save "airfoil_windtunnel_tetrahedra.brep";

/* After creating the geometry we delete and start by loading
 * in the brep file. This is necessary because Gmsh gives the geometries its
 * own indices and after reloading the brep file it uses the brep
 * numeration. */
Delete All;
/* We re-open our brep file. */
Merge "airfoil_windtunnel_tetrahedra.brep";

/* Now, we can create the three-dimensional mesh. */
Mesh.MeshSizeFromCurvature = 4;
Mesh 3;

/* Lastly, we can save the mesh. Note, that we are using msh version 4.X
 * and the parametric option. */
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "airfoil_windtunnel_tetrahedra.msh";
