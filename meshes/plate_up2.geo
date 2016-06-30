
// Geometric parameters
H = 1;  //  semiheight of plate
L = 10; //  width of plate
a = 2;   // crack
b = 6;   // crack

// Discretization parameters
lc1 = .2; // element size at the border
lc2 = .2; // element size at the crack

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {0.0,H,0.0,lc1};
Point(5) = {a,0.0,0.0,lc2};
Point(6) = {b,0.0,0.0,lc2};

Line(1) = {1,5};
Line(2) = {5,6};
Line(3) = {6,2};
Line(4) = {2,3};
Line(5) = {3,4};
Line(6) = {4,1};

Line Loop(11) = {1,2,3,4,5,6};
Plane Surface(1) = {11};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(9) = {1,2,3};
Physical Surface(7) = {1};
