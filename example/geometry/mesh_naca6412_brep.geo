// Gmsh project created on Mon Dec 20 14:25:25 2021
SetFactory("OpenCASCADE");
Merge "naca6412_2.brep";
Transfinite Curve {20, 19, 30, 31, 10, 9, 4, 3, 25, 26, 36, 35} = 6 Using Progression 1;
//+
Transfinite Curve {23, 22, 28, 27, 34, 33} = 4 Using Progression 1;
//+
Transfinite Curve {14, 15, 37, 38, 39, 40, 17, 18} = 9 Using Progression 1;
//+
Transfinite Curve {6, 7, 12, 11} = 6 Using Progression 1;
//+
Transfinite Curve {1, 8, 21, 16, 13, 24, 29, 32, 2, 5} = 2 Using Progression 1;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
Mesh 3;
