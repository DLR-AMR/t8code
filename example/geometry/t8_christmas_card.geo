SetFactory("OpenCASCADE");
Point(1) = {29, 0, 0, 1.0};
Point(2) = {29, 136, 0, 1.0};
Point(3) = {90, 105, 0, 1.0};
Point(4) = {172, 96, 0, 1.0};
Point(5) = {145, 112, 0, 1.0};
Point(6) = {126, 127, 0, 1.0};
Point(7) = {192, 110, 0, 1.0};
Point(8) = {273, 111, 0, 1.0};
Point(9) = {167, 178, 0, 1.0};
Point(10) = {80, 279, 0, 1.0};
Point(11) = {147, 264, 0, 1.0};
Point(12) = {216, 265, 0, 1.0};
Point(13) = {114, 335, 0, 1.0};
Point(14) = {55, 437, 0, 1.0};
Point(15) = {107, 421, 0, 1.0};
Point(16) = {163, 421, 0, 1.0};
Point(17) = {73, 508, 0, 1.0};
Point(18) = {25, 597, 0, 1.0};
Point(19) = {68, 583, 0, 1.0};
Point(20) = {113, 581, 0, 1.0};
Point(21) = {48, 645, 0, 1.0};
Point(22) = {0, 730, 0, 1.0};
Point(23) = {-29, 0, 0, 1.0};
Line(1) = {23, 1};
Line(2) = {1, 2};
Spline(3) = {2:4};
Spline(4) = {4:6};
Spline(5) = {6:8};
Spline(6) = {8:10};
Spline(7) = {10:12};
Spline(8) = {12:14};
Spline(9) = {14:16};
Spline(10) = {16:18};
Spline(11) = {18:20};
Spline(12) = {20:22};
Symmetry {1, 0, 0, 0} {
  Duplicata { Curve{12}; Curve{11}; Curve{10}; Curve{9}; Curve{8}; Curve{7}; Curve{6}; Curve{5}; Curve{4}; Curve{3}; Curve{2};}
}
Coherence;
Curve Loop(1) = {15, 14, 13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 23, 22, 21, 20, 19, 18, 17, 16};
Plane Surface(1) = {1};
Coherence;

Save "t8_christmas_card.brep";
Delete All;

Merge "t8_christmas_card.brep";
Transfinite Curve {:} = 2 Using Progression 1;
Transfinite Curve {16, 14} = 3 Using Progression 1;
Mesh.RecombineAll = 1;
Mesh 2;
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "t8_christmas_card.msh";

