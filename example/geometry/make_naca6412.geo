SetFactory("OpenCASCADE");
// NACA6412 Aerofoil Shape and Spline fit

// We define the points for the spline profile
lc = 0.5;

// trailing edge
Point(1000) = {1.00400,  0.00000, 0, lc};

// dorsal side
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

// leading edge
Point(2000) = {0.00000,  0.00000, 0, lc};

// ventral side
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

// Spline
Spline(100) = {1000:1030, 2000:2028, 1000};

// Points for the flow domain
Point(2032) = {-0.2, -0.5, 0, 1.0};
Point(2033) = {-0.2, 0, 0, 1.0};
Point(2034) = {-0.2, 0.5, 0, 1.0};
Point(2035) = {1.5, 0.5, 0, 1.0};
Point(2036) = {1.5, -0.5, 0, 1.0};
Point(2037) = {1.5, 0, 0, 1.0};
Point(2038) = {1.1, 0.5, 0, 1.0};
Point(2039) = {1.1, -0.5, 0, 1.0};
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
BooleanFragments{ Curve{100}; Delete; }{ Curve{107}; Curve{106}; Curve{108}; }

// Here we define our surfaces and extrude them to get a volume 
Curve Loop(1) = {107, -117, 105, -109};
Plane Surface(1) = {1};
Curve Loop(2) = {110, 111, -103, 105};
Plane Surface(2) = {2};
Curve Loop(3) = {112, 113, -104, 103};
Plane Surface(3) = {3};
Curve Loop(4) = {120, 104, 114, -108};
Curve Loop(5) = {104, 114, -108, 120};
Plane Surface(4) = {5};
Curve Loop(6) = {108, 115, -106, 119};
Plane Surface(5) = {6};
Curve Loop(7) = {116, 107, 118, 106};
Plane Surface(6) = {7};
Extrude {0, 0, 0.1} {
  Surface{1}; Surface{6}; Surface{5}; Surface{4}; Surface{3}; Surface{2}; 
}

// We save the geometry in the brep format
Save "naca6412.brep";

// After creating the geometry we delete everything and start by loading in the brep file.
// This is necessay, because gmsh gives the geometries its own label and after reloading the brep file
// it uses the brep numeration
Delete All;
// We open our brep file
Merge "naca6412.brep";

// After setting our geometries to transfinite we can mesh them
Transfinite Curve {:} = 3 Using Progression 1;
Transfinite Curve {8, 5, 29, 44, 36, 37, 1, 13, 22, 16, 2, 21} = 2 Using Progression 1;
Transfinite Curve {6, 7, 11, 12, 32, 33, 34, 35} = 6 Using Progression 1;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
Mesh 3;

// Now we can save the parametric mesh
Mesh.MshFileVersion = 4.1;
Mesh.SaveParametric = 1;
Save "naca6412.msh";
