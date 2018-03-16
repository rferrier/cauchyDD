
// Geometric parameters
H = 1;  //  height of plate
L = 1;  //  width of plate
a = .4; //  Center of the inclusion
b = .3; //
R = .2; // Radius of the inclusion

// Discretization parameters
lc1 = .05; // element size at the border
lc2 = .05; // element size at the crack tip

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {0.0,H,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Inclusion
Point(5) = {a,b,0.0,lc2};
Point(6) = {a+R,b,0.0,lc2};
Point(7) = {a,b+R,0.0,lc2};
Point(8) = {a-R,b,0.0,lc2};
Point(9) = {a,b-R,0.0,lc2};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5,6,7,8};
Plane Surface(1) = {11,-12};
Plane Surface(2) = {12};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Surface(1) = {1};
Physical Surface(2) = {2};
