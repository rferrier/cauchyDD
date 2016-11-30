// Geometric parameters
H = 20;  //  height of plate
L = 3;  //  width of plate
eps = L/10; // Localization parameter

// Discretization parameters
lc1 = .1; // element size at the border
lc3 = .5; // element size in the domain

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};

// Middle line
Point(5) = {0.0,H/6,0.0,lc1};
Point(6) = {L,H/6,0.0,lc1};

// Refining points
// Point(25) = {0.0,H/6-eps,0.0,lc3};
// Point(26) = {L,H/6-eps,0.0,lc3};

Line(1) = {1,2};
Line(7) = {2,6};
Line(5) = {6,5};
Line(6) = {5,1};

// Refining lines
//Line(16) = {5,25};
//Line(17) = {26,6};
//Line(25) = {25,26};

Line Loop(11) = {1,7,5,6};
Plane Surface(1) = {11};

// Include the refining line
// Line{25} In Surface{1};

Physical Line(1) = {1};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Surface(7) = {1};
