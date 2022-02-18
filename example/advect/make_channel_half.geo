SetFactory("OpenCASCADE");

l = 2;
r1 = 0.5;
r2 = 0.6;
h = 3;
b = 0.2;

Point(1) = {-l - r2, 0, 0};
Point(2) = {l + r2, 0, 0};
Point(3) = {l + r2, h, 0};
Point(4) = {-l - r2, h, 0};

Point(5) = {0, 0, 0};

Point(6) = {Sin(-Pi/2)*r2, Cos(-Pi/2)*r2, 0};
Point(7) = {Sin(-Pi/2)*r1, Cos(-Pi/2)*r1, 0};

Point(8) = {Sin(-Pi/4)*r2, Cos(-Pi/4)*r2, 0};
Point(9) = {Sin(-Pi/4)*r1, Cos(-Pi/4)*r1, 0};

Point(10) = {Sin(Pi/4)*r2, Cos(Pi/4)*r2, 0};
Point(11) = {Sin(Pi/4)*r1, Cos(Pi/4)*r1, 0};

Point(12) = {Sin(Pi/2)*r2, Cos(Pi/2)*r2, 0};
Point(13) = {Sin(Pi/2)*r1, Cos(Pi/2)*r1, 0};

Point(14) = {-l - r2, Cos(-Pi/4)*r2, 0};
Point(15) = {l + r2, Cos(Pi/4)*r2, 0};
Point(16) = {Sin(-Pi/4)*r2, h, 0};
Point(17) = {Sin(Pi/4)*r2, h, 0};

Line(1) = {4, 14};
Line(2) = {14, 1};
Line(3) = {1, 6};
Line(4) = {6, 7};
Line(5) = {13, 12};
Line(6) = {12, 2};
Line(7) = {2, 15};
Line(8) = {15, 3};
Line(9) = {3, 17};
Line(10) = {17, 16};
Line(11) = {16, 4};
Line(12) = {8, 16};
Line(13) = {8, 14};
Line(14) = {10, 17};
Line(15) = {10, 15};
Line(16) = {9, 8};
Line(17) = {11, 10};
Circle(18) = {7, 5, 9};
Circle(19) = {9, 5, 11};
Circle(20) = {11, 5, 13};
Circle(21) = {6, 5, 8};
Circle(22) = {8, 5, 10};
Circle(23) = {10, 5, 12};
Curve Loop(1) = {1, -13, 12, 11};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 21, 13};
Plane Surface(2) = {2};
Curve Loop(3) = {16, -21, 4, 18};
Plane Surface(3) = {3};
Curve Loop(4) = {22, -17, -19, 16};
Plane Surface(4) = {4};
Curve Loop(5) = {14, 10, -12, 22};
Plane Surface(5) = {5};
Curve Loop(6) = {9, -14, 15, 8};
Plane Surface(6) = {6};
Curve Loop(7) = {23, 6, 7, -15};
Plane Surface(7) = {7};
Curve Loop(8) = {17, 23, -5, -20};
Plane Surface(8) = {8};
Extrude {0, 0, b} {
  Surface{1}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{4}; Surface{3}; Surface{2}; 
}

Save "channel_half.brep";
Delete All;
Merge "channel_half.brep";

l_elems = 6;
r2_elems = 6;
h1_elems = 2;
h2_elems = 6;
b_elems = 1;

Transfinite Curve {59, 51, 56, 43, 30, 35, 36, 31, 44, 52, 60, 57} = h1_elems + 1 Using Progression 1;
Transfinite Curve {54, 48, 38, 41, 42, 39, 49, 55} = r2_elems + 1 Using Progression 1;
Transfinite Curve {46, 19, 47, 20, 18, 17} = 2*h1_elems + 1 Using Progression 1;
Transfinite Curve {12, 11, 7, 62, 6, 61} = l_elems + 1 Using Progression 1;
Transfinite Curve {23, 22, 26, 34, 25, 33} = l_elems + 1 Using Progression 1;
Transfinite Curve {3, 4, 10, 9, 16, 15, 28, 27} = h2_elems + 1 Using Progression 1;
Transfinite Curve {8, 1, 2, 14, 58, 5, 50, 21, 45, 53, 13, 13, 13, 13, 37, 40, 29, 24, 32} = b_elems + 1 Using Progression 1;

// After setting our geometries to transfinite we can mesh them
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
Mesh 3;

Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "channel_half.msh";
