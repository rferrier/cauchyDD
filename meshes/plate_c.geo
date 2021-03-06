
// Geometric parameters
H = .2;  //  height of plate
L = 1;  //  width of plate

// Coords of the tops of the crack
//x5 = L/5; y5 = 2*H/3;
///x5 = 2*L/5; y5 = 2*H/3;
x5 = L/5; y5 = 2*H/3;
//x6 = L/5; y6 = 3*H/4;
///x6 = 3*L/4; y6 = 3*H/4;
x6 = L/2; y6 = 3*H/4;

//x5 = L/2; y5 = H/2;
//x6 = 3*L/4; y6 = 5*H/8;
Ri = 100*H;  // Hack : radius of the crack
xm = (x5+x6)/2; ym = (y5+y6)/2;
x7 = xm + Ri; y7 = ym - Ri*(x6-x5)/(y6-y5);
x8 = xm - Ri; y8 = ym + Ri*(x6-x5)/(y6-y5);

// Discretization parameters
lc1 = .02; // element size at the border
lc2 = .02; // element size at the crack tip

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

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Circle(5) = {5,7,6};
Circle(6) = {6,8,5};

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5,6};
Plane Surface(1) = {11,12};

// Include the crack in the mesh
//Point{5} In Surface{1};
//Point{6} In Surface{1};
//Line{5} In Surface{1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Surface(7) = {1};

// Split nodes (this is still not working)
//Plugin(Crack).Dimension = 1 ; 
//Plugin(Crack).PhysicalGroup = (5) ; 
//Plugin(Crack).OpenBoundaryPhysicalGroup = (6) ; 
//Plugin(Crack).Run ;
