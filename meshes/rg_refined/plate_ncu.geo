
// Geometric parameters
H   = 1;  //  height of plate
L   = 1;  //  width of plate
eps = L/15; // refined boundary

// Coords of the tops of the crack

x5 = 2*L/5; y5 = 2*H/3;
x6 = 3*L/4; y6 = H/4;
//x5 = L/5; y5 = 2*H/3;
//x6 = L/2; y6 = 3*H/4;

Ri = 100*H;  // Hack : radius of the crack
xm = (x5+x6)/2; ym = (y5+y6)/2;
x7 = xm + Ri; y7 = ym - Ri*(x6-x5)/(y6-y5);
x8 = xm - Ri; y8 = ym + Ri*(x6-x5)/(y6-y5);

// Discretization parameters
lc1 = .1;  // element size at the interior
lc2 = .1;  // element size at the crack tip
lc3 = .1; // element size at the boundary

// Domain construction
Point(1) = {0.0,0.0,0.0,lc3};
Point(2) = {L,0.0,0.0,lc3};
Point(3) = {L,H,0.0,lc3};
Point(4) = {0.0,H,0.0,lc3};

Point(11) = {eps,eps,0.0,lc1};
Point(12) = {L-eps,eps,0.0,lc1};
Point(13) = {L-eps,H-eps,0.0,lc1};
Point(14) = {eps,H-eps,0.0,lc1};

// Crack tops
Point(5) = {x5, y5, 0.0, lc2};
Point(6) = {x6, y6, 0.0, lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};

Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,11};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

Line{5} In Surface{1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(7) = {1};

