// Geometric parameters
H = 20;  //  height of plate
L = 5;  //  width of plate

// Discretization parameters
lc1 = .5; // element size at the border

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};

// Middle line
Point(5) = {0.0,H/6,0.0,lc1};
Point(6) = {L,H/6,0.0,lc1};

Line(1) = {1,2};
Line(7) = {2,6};
Line(5) = {6,5};
Line(6) = {5,1};

Line Loop(11) = {1,7,5,6};
Plane Surface(1) = {11};

Physical Line(1) = {1};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Surface(7) = {1};
