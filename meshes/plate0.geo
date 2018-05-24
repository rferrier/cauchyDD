
// Geometric parameters
H = 1;  //  height of plate
L = 5;  //  width of plate

// Discretization parameters
lc1 = .1; // element size at the border
lc2 = .1; // element size at the crack tip

// Domain construction
Point(1) = {-L,-H,0.0,lc1};
Point(2) = {L,-H,0.0,lc1};
Point(3) = {L,H,0.0,lc2};
Point(4) = {-L,H,0.0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Surface(5) = {1};
