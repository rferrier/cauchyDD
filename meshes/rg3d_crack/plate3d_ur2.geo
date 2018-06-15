// Geometric parameters
H   = 1;     //  height of plate
L   = 1;     //  width of plate
E   = 1;
// Discretization parameters

lc1 = .05;     // element size at the border
lc2 = .1;

// Domain construction
Point(1) = {0.0,0.0,0.0,lc2};
Point(2) = {L,0.0,0.0,lc2};
Point(3) = {L,E,0.0,lc2};
Point(4) = {0.0,E,0.0,lc2};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,E,H,lc1};
Point(8) = {0.0,E,H,lc1};

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

// Twin domain construction
Point(11) = {0.0,0.0,0.0,lc2};
Point(12) = {L,0.0,0.0,lc2};
Point(13) = {L,E,0.0,lc2};
Point(14) = {0.0,E,0.0,lc2};
Point(15) = {0.0,0.0,H,lc2};
Point(16) = {L,0.0,H,lc2};
Point(17) = {L,E,H,lc2};
Point(18) = {0.0,E,H,lc2};

Line(101) = {11,12};
Line(102) = {12,13};
Line(103) = {13,14};
Line(104) = {14,11};

Line(105) = {16,15};
Line(106) = {17,16};
Line(107) = {18,17};
Line(108) = {15,18};

Line(109) = {11,15};
Line(110) = {12,16};
Line(111) = {13,17};
Line(112) = {14,18};

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

Line Loop(1101) = {101,102,103,104};
Plane Surface(11) = {1101};
Line Loop(1102) = {105,106,107,108};
Plane Surface(12) = {1102};
Line Loop(1103) = {101,110,105,-109};
Plane Surface(13) = {1103};
Line Loop(1104) = {102,111,106,-110};
Plane Surface(14) = {1104};
Line Loop(1105) = {103,112,107,-111};
Plane Surface(15) = {1105};
Line Loop(1106) = {104,109,108,-112};
Plane Surface(16) = {1106};

Surface Loop(1001) = {1,2,3,4,5,6};
Volume(1) = {1001};
Surface Loop(2001) = {11,12,13,14,15,16};
Volume(2) = {2001};

Physical Surface(1) = {1};
Physical Surface(2) = {12};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};

Physical Volume(1) = {1};
Physical Volume(2) = {2};
