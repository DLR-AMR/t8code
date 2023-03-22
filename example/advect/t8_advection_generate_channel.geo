//  This file is part of t8code.
//  t8code is a C library to manage a collection (a forest) of multiple
//  connected adaptive space-trees of general element types in parallel.
//
//  Copyright (C) 2015 the developers
//
//  t8code is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  t8code is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with t8code; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

// This file will create a brep and msh file of a 
// cuboidal flow channel with a cylindrical hole in the middle
// when opened with Gmsh.

// Use the OpenCASCADE geometry kernel to generate a brep file later on
SetFactory("OpenCASCADE");

// Set the dimensions of the channel
r1 = 0.5; // Radius of the hole
r2 = 0.6; // Radius of the ring around the hole
h = 2; // height/2 of the channel
l = 3; // length/2 of the channel
b = 3; // width of the channel

// Definition all the points needed for construction
Point(1) = {0, 0, -b/2};

Point(2) = {Sin(-Pi+Pi/4)*r2, Cos(-Pi+Pi/4)*r2, -b/2};
Point(3) = {Sin(-Pi+Pi/4)*r1, Cos(-Pi+Pi/4)*r1, -b/2};

Point(4) = {Sin(-Pi/4)*r2, Cos(-Pi/4)*r2, -b/2};
Point(5) = {Sin(-Pi/4)*r1, Cos(-Pi/4)*r1, -b/2};

Point(6) = {Sin(Pi/4)*r2, Cos(Pi/4)*r2, -b/2};
Point(7) = {Sin(Pi/4)*r1, Cos(Pi/4)*r1, -b/2};

Point(8) = {Sin(Pi-Pi/4)*r2, Cos(Pi-Pi/4)*r2, -b/2};
Point(9) = {Sin(Pi-Pi/4)*r1, Cos(Pi-Pi/4)*r1, -b/2};

Point(10) = {-l, h, -b/2};
Point(11) = {-l, -h, -b/2};
Point(12) = {l, h, -b/2};
Point(13) = {l, -h, -b/2};

Point(14) = {-l, Cos(-Pi+Pi/4)*r2, -b/2};
Point(15) = {-l, Cos(-Pi/4)*r2, -b/2};
Point(16) = {-l, Cos(Pi/4)*r2, -b/2};
Point(17) = {-l, Cos(Pi-Pi/4)*r2, -b/2};

Point(18) = {l, Cos(-Pi+Pi/4)*r2, -b/2};
Point(19) = {l, Cos(-Pi/4)*r2, -b/2};
Point(20) = {l, Cos(Pi/4)*r2, -b/2};
Point(21) = {l, Cos(Pi-Pi/4)*r2, -b/2};

Point(22) = {Sin(-Pi+Pi/4)*r2, h, -b/2};
Point(23) = {Sin(-Pi/4)*r2, h, -b/2};
Point(24) = {Sin(Pi/4)*r2, h, -b/2};
Point(25) = {Sin(Pi-Pi/4)*r2, -b/2};

Point(26) = {Sin(-Pi+Pi/4)*r2, -h, -b/2};
Point(27) = {Sin(-Pi/4)*r2, -h, -b/2};
Point(28) = {Sin(Pi/4)*r2, -h, -b/2};
Point(29) = {Sin(Pi-Pi/4)*r2, -h, -b/2};

// Definition all the one dimensional geometries
Circle(1) = {6, 1, 8};
Circle(2) = {7, 1, 9};
Circle(3) = {8, 1, 2};
Circle(4) = {9, 1, 3};
Circle(5) = {2, 1, 4};
Circle(6) = {3, 1, 5};
Circle(7) = {4, 1, 6};
Circle(8) = {5, 1, 7};
Line(9) = {5, 4};
Line(10) = {7, 6};
Line(11) = {9, 8};
Line(12) = {3, 2};
Line(13) = {10, 22};
Line(14) = {22, 24};
Line(15) = {24, 12};
Line(16) = {12, 19};
Line(17) = {19, 18};
Line(18) = {18, 13};
Line(19) = {13, 28};
Line(20) = {28, 26};
Line(21) = {26, 11};
Line(22) = {11, 14};
Line(23) = {14, 15};
Line(24) = {15, 10};
Line(25) = {6, 24};
Line(26) = {6, 19};
Line(27) = {18, 8};
Line(28) = {8, 28};
Line(29) = {2, 26};
Line(30) = {2, 14};
Line(31) = {15, 4};
Line(32) = {4, 22};

// Definition of the surfaces
Curve Loop(1) = {31, -5, 30, 23};
Curve Loop(2) = {31, -5, 30, 23};
Plane Surface(1) = {2};
Curve Loop(3) = {24, 13, -32, -31};
Plane Surface(2) = {3};
Curve Loop(4) = {25, -14, -32, 7};
Plane Surface(3) = {4};
Curve Loop(5) = {26, -16, -15, -25};
Plane Surface(4) = {5};
Curve Loop(6) = {27, -1, 26, 17};
Curve Loop(7) = {26, 17, 27, -1};
Plane Surface(5) = {7};
Curve Loop(8) = {18, 19, -28, -27};
Plane Surface(6) = {8};
Curve Loop(9) = {3, 29, -20, -28};
Plane Surface(7) = {9};
Curve Loop(10) = {29, 21, 22, -30};
Plane Surface(8) = {10};
Curve Loop(11) = {6, 9, -5, -12};
Curve Loop(12) = {6, 9, -5, -12};
Plane Surface(9) = {12};
Curve Loop(13) = {8, 10, -7, -9};
Plane Surface(10) = {13};
Curve Loop(14) = {2, 11, -1, -10};
Curve Loop(15) = {2, 11, -1, -10};
Plane Surface(11) = {15};
Curve Loop(16) = {3, -12, -4, 11};
Plane Surface(12) = {16};

// Extrusion of the two dimensional surfaces
Extrude {0, 0, b} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{10}; Surface{11}; Surface{12}; 
}
// Make sure, that no geometries are defined multiple times
Coherence;

// Save the channel geometry and reopen it. 
// Gmsh has its own numbering of the geometries during the construction, 
// but we need the numbering inside the generated brep file. 
// Therefore, we delete everything and open the geometry inside the brep file.
Save "channel.brep";
Delete All;
Merge "channel.brep";

// Definition of the amount of elements in each direction
l_elems = 5; // Amount of elements/2 along the length of the channel 
h_elems = 3; // Amount of elements/2 along the height of the channel
r2_elems = 1; // Amount of "boundary layer" elements around the hole
circle_elems = 2; // Amount of elements/4 around the hole
b_elems = 5; // Amount of elements along the width of the channel

// After setting our geometries to transfinite we can mesh them
Transfinite Curve {61, 62, 9, 10, 3, 4, 17, 18, 35, 36, 31, 30, 42, 41, 50, 49} = l_elems + 1 Using Progression 1;
Transfinite Curve {63, 64, 56, 57, 51, 52, 46, 47, 33, 34, 23, 24, 19, 20, 14, 15} = h_elems + 1 Using Progression 1;
Transfinite Curve {69, 70, 72, 71, 81, 82, 77, 76} = r2_elems + 1 Using Progression 1;
Transfinite Curve {27, 28, 74, 75, 67, 67, 68, 6, 7, 53, 54, 83, 84, 79, 80, 43, 43, 44, 11, 12, 26, 25, 59, 58, 38, 39} = circle_elems + 1 Using Progression 1;
Transfinite Curve {16, 13, 1, 8, 60, 22, 2, 66, 5, 65, 73, 21, 40, 78, 55, 48, 32, 29, 37, 45} = b_elems + 1 Using Progression 1;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
Mesh 3;
Coherence Mesh;

// Save the mesh
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "channel.msh";
