
// Geometric parameters
H   = 100;  //  height of plate
L   = 100;  //  width of plate
U   = 20;   // Loading place
eps = L/30; // refined boundary
ep2 = L/2000;

// Discretization parameters
lc1 = 5;  // element size at the interior
lc2 = 2.5;  // element size at the crack tip
lc3 = 5; // element size at the boundary

// Coords of the tops of the crack
x5 = L/2;   y5 = H/2;
x6 = L/2+ep2/2;   y6 = H-U;
x7 = L/2-ep2/2;   y7 = H-U;

// Domain construction
Point(1) = {0.0,0.0,0.0,lc3};
Point(2) = {L,0.0,0.0,lc3};
Point(3) = {L,H,0.0,lc3};
Point(4) = {0.0,H,0.0,lc3};

Point(8) = {L/2+U/2,H,0.0,lc3};
Point(9) = {L/2-U/2,H,0.0,lc3};
Point(10) = {L/2+U/2,H-U,0.0,lc3};
Point(11) = {L/2-U/2,H-U,0.0,lc3};

//Point(11) = {eps,eps,0.0,lc1};
//Point(12) = {L-eps,eps,0.0,lc1};
//Point(13) = {L-eps,H-eps,0.0,lc1};
//Point(14) = {eps,H-eps,0.0,lc1};

// Crack tops
Point(5) = {x5, y5, 0.0, lc2};
Point(6) = {x6, y6, 0.0, lc2};
Point(7) = {x7, y7, 0.0, lc2};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,8};
Line(4)  = {8,10};
Line(5)  = {10,6};
Line(6)  = {6,5};
Line(7)  = {5,7};
Line(8)  = {7,11};
Line(9)  = {11,9};
Line(10) = {9,4};
Line(11) = {4,1};

//Line(11) = {11,12};
//Line(12) = {12,13};
//Line(13) = {13,14};
//Line(14) = {14,11};

Line Loop(11) = {1,2,3,4,5,6,7,8,9,10,11};
Plane Surface(1) = {11};

//Line{11} In Surface{1};
//Line{12} In Surface{1};
//Line{13} In Surface{1};
//Line{14} In Surface{1};

Physical Line(1)  = {1};
Physical Line(2)  = {2};
Physical Line(3)  = {3};
Physical Line(4)  = {4}; 
Physical Line(5)  = {5};
Physical Line(6)  = {6};
Physical Line(7)  = {7};
Physical Line(8)  = {8};
Physical Line(9)  = {9};
Physical Line(10) = {10};
Physical Line(11) = {11};

Physical Surface(1) = {1};

