
// Geometric parameters
H = 7;  //  semiheight of plate
L = 5;  //  semiwidth of plate
a = 1; // semilength of crack (center of crack is at x1=x2=0)

// Discretization parameters
lc1 = .5; // element size at the border
lc2 = .1; // element size at the top

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc2};
Point(4) = {0.0,H,0.0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4}; 
Physical Surface(5) = {1};
