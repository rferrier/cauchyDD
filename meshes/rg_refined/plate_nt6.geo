
// Geometric parameters
H   = 15;  //  height of plate
L   = 10;  //  width of plate
eps = L/30; // refined boundary

// Discretization parameters
lc1 = .5; // element size at the interior
lc2 = .1;

// Domain construction
Point(1) = {0.0,0.0,0.0,lc2};
Point(2) = {L,0.0,0.0,lc2};
Point(3) = {L,H,0.0,lc2};
Point(4) = {0.0,H,0.0,lc2};

Point(11) = {eps,eps,0.0,lc1};
Point(12) = {L-eps,eps,0.0,lc1};
Point(13) = {L-eps,H-eps,0.0,lc1};
Point(14) = {eps,H-eps,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,11};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

Line{11} In Surface{1};
Line{12} In Surface{1};
Line{13} In Surface{1};
Line{14} In Surface{1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Surface(5) = {1};
