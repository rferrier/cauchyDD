// Geometric parameters
H   = 20;  //  height of plate
L   = 5;   //  width of plate
eps = L/10; // Localization parameter

// Coords of the tops of the crack
//x5 = L/5; y5 = 2*H/3;
x7 = 4*L/10; y7 = 2*H/4;
x8 = 5*L/10; y8 = 3*H/4;
//x6 = 3*L/4; y6 = 3*H/4;
//x5 = L/2; y5 = H/2;
//x6 = 3*L/4; y6 = 5*H/8;
Ri = 100*H;  // Hack : radius of the crack
xm = (x7+x8)/2; ym = (y7+y8)/2;
x9 = xm + Ri; y9 = ym - Ri*(x8-x7)/(y8-y7);
x10 = xm - Ri; y10 = ym + Ri*(x8-x7)/(y8-y7);

// Discretization parameters
lc1 = .05; // element size at the border
lc2 = .5; // element size at the crack tip
lc3 = .5; // element size in the domain

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {0.0,H,0.0,lc1};

// Middle line
Point(5) = {0.0,H/6,0.0,lc1};
Point(6) = {L,H/6,0.0,lc1};

// Crack tops
Point(7) = {x7, y7, 0.0, lc2};
Point(8) = {x8, y8, 0.0, lc2};
Point(9) = {x9, y9, 0.0, lc2};
Point(10) = {x10, y10, 0.0, lc2};

// Refining points
Point(13) = {L-eps,H-eps,0.0,lc3};
Point(14) = {eps,H-eps,0.0,lc3};
Point(15) = {eps,H/6+eps,0.0,lc3};
Point(16) = {L-eps,H/6+eps,0.0,lc3};

Point(21) = {eps,eps,0.0,lc3};
Point(22) = {L-eps,eps,0.0,lc3};
Point(26) = {L-eps,H/6-eps,0.0,lc3};
Point(25) = {eps,H/6-eps,0.0,lc3};

Line(1) = {1,2};
Line(7) = {2,6};
Line(2) = {6,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {6,5};
Line(6) = {5,1};
Circle(8) = {7,9,8};
Circle(9) = {8,10,7};

// Refining lines
Line(12) = {16,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(21) = {21,22};
Line(27) = {22,26};
Line(25) = {26,25};
Line(26) = {25,21};

Line Loop(11) = {1,7,2,3,4,6};
Line Loop(12) = {8,9};
Plane Surface(1) = {11,12};

// Include the middle line in the mesh
Line{5} In Surface{1};

// Include the refining lines
Line{12} In Surface{1};
Line{13} In Surface{1};
Line{14} In Surface{1};
Line{15} In Surface{1};
Line{21} In Surface{1};
Line{27} In Surface{1};
Line{25} In Surface{1};
Line{26} In Surface{1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Surface(7) = {1};
