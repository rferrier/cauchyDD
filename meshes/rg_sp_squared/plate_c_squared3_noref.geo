// Geometric parameters
H   = 1;  //  height of plate
L   = 1;  //  width of plate
eps = H/30; // refining width

// Coords of the tips of the crack
x5 = L/6;   y5 = H/9;//
x6 = 3*L/7;  y6 = H/7; //x6 = L/2; y6 = H/5;//

Ri = 100*H;  // Hack : radius of the crack
xm = (x5+x6)/2; ym = (y5+y6)/2;
x7 = xm + Ri; y7 = ym - Ri*(x6-x5)/(y6-y5);
x8 = xm - Ri; y8 = ym + Ri*(x6-x5)/(y6-y5);

// Discretization parameters
lc1 = .002;  // element size at the boundary
lc2 = .025;  // element size at the crack tip
lc3 = .05;  // element size at the interior

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {0.0,H,0.0,lc1};

// Crack tops
Point(5) = {x5, y5, 0.0, lc2};
Point(6) = {x6, y6, 0.0, lc2};
Point(7) = {x7, y7, 0.0, lc2};
Point(8) = {x8, y8, 0.0, lc2};

// Middle line
Point(9)  = {0.0,H/2,0.0,lc1};
Point(10) = {L,H/2,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,10};
Line(3) = {10,9};
Line(4) = {9,1};
Circle(5) = {5,7,6};
Circle(6) = {6,8,5};

Line(7) = {10,3};
Line(8) = {3,4};
Line(9) = {4,9};

// Refining lines
Point(11) = {eps,eps,0.0,lc3};
Point(12) = {L-eps,eps,0.0,lc3};
Point(13) = {L-eps,H/2-eps,0.0,lc3};
Point(14) = {eps,H/2-eps,0.0,lc3};

Point(15) = {eps,H/2+eps,0.0,lc3};
Point(16) = {L-eps,H/2+eps,0.0,lc3};
Point(17) = {L-eps,H-eps,0.0,lc3};
Point(18) = {eps,H-eps,0.0,lc3};

Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,11};

Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,15};

//

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5,6};
Line Loop(13) = {-3,7,8,9};
Plane Surface(1) = {11,12};
Plane Surface(2) = {13};

Line{11} In Surface{1};
Line{12} In Surface{1};
Line{13} In Surface{1};
Line{14} In Surface{1};

Line{15} In Surface{2};
Line{16} In Surface{2};
Line{17} In Surface{2};
Line{18} In Surface{2};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7}; 
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Surface(1) = {1};
Physical Surface(2) = {2};

