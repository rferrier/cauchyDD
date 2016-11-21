// Geometric parameters
H = 20;  //  height of plate
L = 5;  //  width of plate

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
lc1 = .5; // element size at the border
lc2 = .1; // element size at the crack tip

// Domain construction
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

Line(2) = {6,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Circle(8) = {7,9,8};
Circle(9) = {8,10,7};

Line Loop(11) = {2,3,4,5};
Line Loop(12) = {8,9};
Plane Surface(1) = {11,12};

// Include the middle line in the mesh
Line{5} In Surface{1};

Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Line(5) = {5};
Physical Surface(7) = {1};
