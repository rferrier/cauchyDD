// Geometric parameters
H   = 1;     //  height of plate
L   = 1;     //  width of plate
E   = 1;
// Discretization parameters

lc1 = .025;     // element size at the border
lc2 = .1;

// Domain construction
Point(1)  = {0.0,0.0,0.0,lc2};
Point(2)  = {L,0.0,0.0,lc2};
Point(3)  = {L,E,0.0,lc2};
Point(4)  = {0.0,E,0.0,lc2};
Point(5)  = {0.0,0.0,H,lc2};
Point(6)  = {L,0.0,H,lc2};
Point(7)  = {L,E,H,lc2};
Point(8)  = {0.0,E,H,lc2};
Point(9)  = {0.0,0.0,H-lc1,lc1};
Point(10) = {L,0.0,H-lc1,lc1};
Point(11) = {L,E,H-lc1,lc1};
Point(12) = {0.0,E,H-lc1,lc1};
//Point(91)  = {lc1,0.0,H-lc1,lc1};
//Point(92)  = {0.0,lc1,H-lc1,lc1};
//Point(101) = {L-lc1,0.0,H-lc1,lc1};
//Point(102) = {L,lc1,H-lc1,lc1};
//Point(111) = {L-lc1,E,H-lc1,lc1};
//Point(112) = {L,E-lc1,H-lc1,lc1};
//Point(121) = {lc1,E,H-lc1,lc1};
//Point(122) = {0.0,E-lc1,H-lc1,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {6,5};
Line(6) = {7,6};
Line(7) = {8,7};
Line(8) = {5,8};

Line(91) = {1,9};
Line(101) = {2,10};
Line(111) = {3,11};
Line(121) = {4,12};
Line(92) = {9,5};
Line(102) = {10,6};
Line(112) = {11,7};
Line(122) = {12,8};

Line(13) = {9,10};
Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {12,9};

Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Plane Surface(2) = {102};
Line Loop(103) = {1,101,102,5,-92,-91};
Plane Surface(3) = {103};
Line Loop(104) = {2,111,112,6,-102,-101};
Plane Surface(4) = {104};
Line Loop(105) = {3,121,122,7,-112,-111};
Plane Surface(5) = {105};
Line Loop(106) = {4,91,92,8,-122,-121};
Plane Surface(6) = {106};

Line{13} In Surface{3};
Line{14} In Surface{4};
Line{15} In Surface{5};
Line{16} In Surface{6};

Surface Loop(1001) = {1,2,3,4,5,6};
Volume(1) = {1001};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};

Physical Volume(1) = {1};
