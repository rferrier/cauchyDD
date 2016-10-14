
// Geometric parameters
H = 10;  //  height of plate
L = 5;  //  width of plate

// Discretization parameters
lc1 = 2; // element size at the border
lc2 = 2; // element size at the crack tip

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,L,0.0,lc2};
Point(4) = {0.0,L,0.0,lc2};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,L,H,lc2};
Point(8) = {0.0,L,H,lc2};

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

Surface Loop(1001) = {1,2,3,4,5,6};
Volume(7) = {1001};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};

Physical Volume(7) = {7};
