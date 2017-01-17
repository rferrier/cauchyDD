
// Geometric parameters
H = 2;  //  height of plate
L = 2;  //  width of plate
R = 1;  //  Radius of the hole

// Discretization parameters
lc1 = .1; // element size
lc2 = .1;

// Domain construction
Point(1) = {-L,-H,0.0,lc1};
Point(2) = {L,-H,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {-L,H,0.0,lc1};
// Hole
Point(5) = {0.0,0.0,0.0,lc1};
Point(6) = {R,0.0,0.0,lc1};
Point(7) = {0.0,R,0.0,lc2};
Point(8) = {-R,0.0,0.0,lc1};
Point(9) = {0.0,-R,0.0,lc1};

//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {8,7,6,5};
Plane Surface(1) = {11,12};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Line(5) = {8,7,6,5};
Physical Line(6) = {1,2,4};
Physical Surface(5) = {1};
