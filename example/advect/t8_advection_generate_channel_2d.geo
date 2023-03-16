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
r = 0.5; // Radius of the hole
h = 2; // height/2 of the channel
l = 3; // length/2 of the channel

// Definition all the points needed for construction
Rectangle(1) = {-l, -h, 0, 2*l, 2*h, 0};
Circle(5) = {0, 0, 0, r, 0, 2*Pi};
Curve Loop(2) = {5};
Plane Surface(2) = {2};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

// Save the channel geometry and reopen it. 
// Gmsh has its own numbering of the geometries during the construction, 
// but we need the numbering inside the generated brep file. 
// Therefore, we delete everything and open the geometry inside the brep file.
Save "channel_2d.brep";
Delete All;
Merge "channel_2d.brep";

/* We want a full-quad mesh and therefore have to use the right meshing algorithms.
 * Here we use the Automatic 2D meshing algorithm and the 
 * blossom full-quad recombination algorithm. Feel free to experiment here.
 * The following algorithms are available with Gmsh 4.11.0:
 *
 * 2D mesh algorithms:             1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 
 *                                 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 
 *                                 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 
 *                                11: Quasi-structured Quad
 * Mesh recombination algorithms   0: simple, 1: blossom, 
 *                                 2: simple full-quad, 3: blossom full-quad
 *
 * For other Gmsh versions check the Gmsh website: 
 * https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options */
 Mesh.Algorithm = 2;
 Mesh.RecombineAll = 1;
 Mesh.RecombinationAlgorithm = 3;
 
 /* Now we can create the two-dimensional mesh. */
 Mesh 2;

// Save the mesh
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "channel_2d.msh";
