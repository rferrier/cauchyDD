// Enable cracks (or not...)
Geometry.AutoCoherence = 0;

// Geometric parameters
pi    = 3.1415926535;  // What I love to make learn a number useful to wises
R     = 2;   // Radius of the crack
a     = 4;      // Center of the crack
b     = 3;      //
da    = pi/15; // angle of the crack
H     = 2;  //  height of plate
L     = 7;  //  width of plate
E     = 10; //E     = 7;
eps   = .4;  //refining width

// Discretization parameters
lc1 = .5; // element size at the crack tip
lc2 = .5; // element size in teh domain
lc3 = .2; // refining size

// Domain construction
Point(1) = {0.0,0.0,0.0,lc3};
Point(2) = {L,0.0,0.0,lc3};
Point(3) = {L,E,0.0,lc3};
Point(4) = {0.0,E,0.0,lc3};
Point(5) = {0.0,0.0,H,lc3};
Point(6) = {L,0.0,H,lc3};
Point(7) = {L,E,H,lc3};
Point(8) = {0.0,E,H,lc3};

// Refining points
Point(101) = {eps,eps,eps,lc2};
Point(102) = {L-eps,eps,eps,lc2};
Point(103) = {L-eps,E-eps,eps,lc2};
Point(104) = {eps,E-eps,eps,lc2};
Point(105) = {eps,eps,H-eps,lc2};
Point(106) = {L-eps,eps,H-eps,lc2};
Point(107) = {L-eps,E-eps,H-eps,lc2};
Point(108) = {eps,E-eps,H-eps,lc2};

// aretes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {6,5};
Line(6) = {7,6};
Line(7) = {8,7};
Line(8) = {5,8};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// refining aretes
Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,101};

Line(105) = {106,105};
Line(106) = {107,106};
Line(107) = {108,107};
Line(108) = {105,108};

Line(109) = {101,105};
Line(110) = {102,106};
Line(111) = {103,107};
Line(112) = {104,108};

// External faces
Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Plane Surface(2) = {102};
Line Loop(103) = {1,10,5,-9};
Plane Surface(3) = {103};
Line Loop(104) = {2,11,6,-10};
Plane Surface(4) = {104};
Line Loop(105) = {3,12,7,-11};
Plane Surface(5) = {105};
Line Loop(106) = {4,9,8,-12};
Plane Surface(6) = {106};

// Refining faces
Line Loop(1101) = {101,102,103,104};
Plane Surface(101) = {1101};
Line Loop(1102) = {105,106,107,108};
Plane Surface(102) = {1102};
Line Loop(1103) = {101,110,105,-109};
Plane Surface(103) = {1103};
Line Loop(1104) = {102,111,106,-110};
Plane Surface(104) = {1104};
Line Loop(1105) = {103,112,107,-111};
Plane Surface(105) = {1105};
Line Loop(1106) = {104,109,108,-112};
Plane Surface(106) = {1106};

// Crack construction
Point(9) = {R*Cos(da)+a,b,R*Sin(da)+H/2,lc1};
Point(10) = {R*Cos(pi/2)+a,R*Sin(pi/2)+b,H/2,lc1};
Point(11) = {R*Cos(pi)*Cos(da)+a,R*Sin(pi)+b,-R*Sin(da)+H/2,lc1};
Point(12) = {R*Cos(3*pi/2)+a,R*Sin(3*pi/2)+b,H/2,lc1};
Point(13) = {a,b,H/2,lc1};

Circle(13) = {9,13,10};
Circle(14) = {10,13,11};
Circle(15) = {11,13,12};
Circle(16) = {12,13,9};

Line Loop(107) = {13,14,15,16};
Plane Surface(7) = {107};

// Physical stuff
Surface Loop(1001) = {1,2,3,4,5,6};
Surface Loop(1007) = {7,-8};
Volume(9) = {1001};

// Crack and refining surfaces
Surface{7}   In Volume{9};
Surface{101} In Volume{9};
Surface{102} In Volume{9};
Surface{103} In Volume{9};
Surface{104} In Volume{9};
Surface{105} In Volume{9};
Surface{106} In Volume{9};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};

Physical Volume(9) = {9};
