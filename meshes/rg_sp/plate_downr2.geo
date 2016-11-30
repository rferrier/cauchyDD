// Geometric parameters
H = 20;  //  height of plate
L = 5;  //  width of plate
eps = L/10; // Localization parameter

// Discretization parameters
lc1 = .05; // element size at the border
lc3 = .5; // element size in the domain

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};

// Middle line
Point(5) = {0.0,H/6,0.0,lc1};
Point(6) = {L,H/6,0.0,lc1};

// Refining points
Point(21) = {eps,eps,0.0,lc3};
Point(22) = {L-eps,eps,0.0,lc3};
Point(26) = {L-eps,H/6-eps,0.0,lc3};
Point(25) = {eps,H/6-eps,0.0,lc3};

Line(1) = {1,2};
Line(7) = {2,6};
Line(5) = {6,5};
Line(6) = {5,1};

// Refining lines
Line(21) = {21,22};
Line(27) = {22,26};
Line(25) = {26,25};
Line(26) = {25,21};

Line Loop(11) = {1,7,5,6};
Plane Surface(1) = {11};

// Include the refining line
Line{21} In Surface{1};
Line{27} In Surface{1};
Line{25} In Surface{1};
Line{26} In Surface{1};

Physical Line(1) = {1};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Surface(7) = {1};
