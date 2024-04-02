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

/* This file will create a brep and msh file of a NACA 6412 airfoil.
 * We can then use these files to generate a curved adaptive mesh in
 * t8code.
*/

/* We select the OpenCASCADE geometry kernel because we use the same kernel
 * in t8code. Furthermore, it has a bigger functionality than the inbuilt
 * GEO kernel. */
SetFactory("OpenCASCADE");

/* Now, we can define the points of our airfoil. The coordinates
 * were taken from https://m-selig.ae.illinois.edu/ads/coord/naca6412.dat */
lc = 0.5;

/* Trailing edge */
Point(1000) = {1.00400,  0.00000, 0, lc};

/* Dorsal side */
Point(1001) = {1.00025,  0.00124, 0, lc};
Point(1002) = {0.99758,  0.00216, 0, lc};
Point(1003) = {0.98961,  0.00490, 0, lc};
Point(1004) = {0.97640,  0.00935, 0, lc};
Point(1005) = {0.95808,  0.01538, 0, lc};
Point(1006) = {0.93481,  0.02278, 0, lc};
Point(1007) = {0.90682,  0.03130, 0, lc};
Point(1008) = {0.87436,  0.04068, 0, lc};
Point(1009) = {0.83777,  0.05062, 0, lc};
Point(1010) = {0.79740,  0.06082, 0, lc};
Point(1011) = {0.75366,  0.07097, 0, lc};
Point(1012) = {0.70702,  0.08079, 0, lc};
Point(1013) = {0.65797,  0.08998, 0, lc};
Point(1014) = {0.60703,  0.09827, 0, lc};
Point(1015) = {0.55477,  0.10543, 0, lc};
Point(1016) = {0.50176,  0.11124, 0, lc};
Point(1017) = {0.44863,  0.11554, 0, lc};
Point(1018) = {0.39587,  0.11817, 0, lc};
Point(1019) = {0.34306,  0.11841, 0, lc};
Point(1020) = {0.29199,  0.11583, 0, lc};
Point(1021) = {0.24336,  0.11060, 0, lc};
Point(1022) = {0.19780,  0.10302, 0, lc};
Point(1023) = {0.15592,  0.09344, 0, lc};
Point(1024) = {0.11825,  0.08231, 0, lc};
Point(1025) = {0.08524,  0.07012, 0, lc};
Point(1026) = {0.05726,  0.05736, 0, lc};
Point(1027) = {0.03460,  0.04452, 0, lc};
Point(1028) = {0.01745,  0.03204, 0, lc};
Point(1029) = {0.00595,  0.02029, 0, lc};
Point(1030) = {0.00014,  0.00955, 0, lc};

/* Leading edge */
Point(2000) = {0.00000,  0.00000, 0, lc};

/* Ventral side */
Point(2001) = {0.00534, -0.00792, 0, lc};
Point(2002) = {0.01590, -0.01383, 0, lc};
Point(2003) = {0.03149, -0.01781, 0, lc};
Point(2004) = {0.05186, -0.01999, 0, lc};
Point(2005) = {0.07672, -0.02054, 0, lc};
Point(2006) = {0.10574, -0.01967, 0, lc};
Point(2007) = {0.13861, -0.01763, 0, lc};
Point(2008) = {0.17495, -0.01470, 0, lc};
Point(2009) = {0.21441, -0.01121, 0, lc};
Point(2010) = {0.25664, -0.00748, 0, lc};
Point(2011) = {0.30127, -0.00384, 0, lc};
Point(2012) = {0.34792, -0.00064, 0, lc};
Point(2013) = {0.39622,  0.00182, 0, lc};
Point(2014) = {0.44685,  0.00370, 0, lc};
Point(2015) = {0.49824,  0.00542, 0, lc};
Point(2016) = {0.54976,  0.00684, 0, lc};
Point(2017) = {0.60088,  0.00786, 0, lc};
Point(2018) = {0.65105,  0.00843, 0, lc};
Point(2019) = {0.69972,  0.00853, 0, lc};
Point(2020) = {0.74634,  0.00819, 0, lc};
Point(2021) = {0.79039,  0.00747, 0, lc};
Point(2022) = {0.83137,  0.00643, 0, lc};
Point(2023) = {0.86878,  0.00520, 0, lc};
Point(2024) = {0.90220,  0.00386, 0, lc};
Point(2025) = {0.93121,  0.00252, 0, lc};
Point(2026) = {0.95546,  0.00129, 0, lc};
Point(2027) = {0.97465,  0.00024, 0, lc};
Point(2028) = {1.00000,  0.00000, 0, lc};

/* Build a spline fit for the points. Starting and ending with the 
 * trailing edge. */
Spline(100) = {1000:1030, 2000:2028, 1000};

/* Definition of the corner points of the flow domain. Note, that we
 * to define more than the four corner points. With these additional
 * points we can divide the domain into quadrilaterals and later on
 * into hexahedra. This is important to structurally mesh the domain. */
Point(2032) = {-0.2, -0.5, 0, 1.0};
Point(2033) = {-0.2, 0, 0, 1.0};
Point(2034) = {-0.2, 0.5, 0, 1.0};
Point(2035) = {1.5, 0.5, 0, 1.0};
Point(2036) = {1.5, -0.5, 0, 1.0};
Point(2037) = {1.5, 0, 0, 1.0};
Point(2038) = {1.1, 0.5, 0, 1.0};
Point(2039) = {1.1, -0.5, 0, 1.0};

/* Build lines through the points. */
Line(103) = {1000, 2037};
Line(104) = {1000, 2039};
Line(105) = {1000, 2038};
Line(106) = {2000, 2033};
Line(107) = {2034, 1021};
Line(108) = {2009, 2032};
Line(109) = {2034, 2038};
Line(110) = {2038, 2035};
Line(111) = {2035, 2037};
Line(112) = {2037, 2036};
Line(113) = {2036, 2039};
Line(114) = {2039, 2032};
Line(115) = {2032, 2033};
Line(116) = {2033, 2034};

/* Split the NACA spline into three sections to get quadrilaterals. */
BooleanFragments{ Curve{100}; Delete; }{ Curve{107}; Curve{106}; Curve{108}; }

/* Now, we define our quadrilateral surfaces. */ 
Curve Loop(1) = {107, -117, 105, -109};
Plane Surface(1) = {1};

Curve Loop(2) = {110, 111, -103, 105};
Plane Surface(2) = {2};

Curve Loop(3) = {112, 113, -104, 103};
Plane Surface(3) = {3};

Curve Loop(4) = {104, 114, -108, 120};
Plane Surface(4) = {4};

Curve Loop(5) = {108, 115, -106, 119};
Plane Surface(5) = {5};

Curve Loop(6) = {116, 107, 118, 106};
Plane Surface(6) = {6};

/* And extrude them to get hexahedral volumes. */
Extrude {0, 0, 0.1} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; 
}

/* Due to the indexing behavior if Gmsh, we have to save the .brep file
 * and reopen it again. */
Save "airfoil_windtunnel_hexahedra.brep";

/* After creating the geometry we delete everything and start by loading
 * in the brep file. This is necessary because Gmsh gives the geometries its
 * own indices and after reloading the brep file it uses the brep 
 * numeration. */
Delete All;
/* We re-open our brep file. */
Merge "airfoil_windtunnel_hexahedra.brep";

/* With the new numeration, we can start to transfinite the curves of the
 * geometry. In Gmsh this is the first step of creating a structured mesh.
 * Feel free to change these values and experiment.
 * First, we transfinite all curves with two points (this means, that Gmsh
 * will create two nodes on these curves). */
Transfinite Curve {:} = 2 Using Progression 1;

/* After that, we transfinite all curves, which should have more than two
 * nodes on them. Here we transfinite all curves orthogonal to the airfoil. */
Transfinite Curve {18, 17, 9, 10, 3, 4, 41, 42, 33, 34, 27, 28, 22, 23} = 
  3 Using Progression 1;

/* Then, we transfinite the curves parallel to the airfoil. Because of the difference
 * in the length of the curves, we transfinite some with three nodes and some with
 * six nodes. */
Transfinite Curve {14, 15, 19, 20, 25, 26, 45, 46, 47, 48, 38, 39, 43, 44} = 
  3 Using Progression 1;
Transfinite Curve {11, 12, 6, 7, 35, 36, 30, 31} = 6 Using Progression 1;

/* Then, we can transfinite the surfaces and volumes. Gmsh derives the
 * necessary amount of nodes automatically from the transfinite curves. */
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};

/* Now we can create the three-dimensional mesh. */
Mesh 3;

/* Lastly, we can save the mesh. Note, that we are using msh version 4.X
 * and the parametric option. */
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "airfoil_windtunnel_hexahedra.msh";
