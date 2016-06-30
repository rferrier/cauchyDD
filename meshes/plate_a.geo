
// Geometric parameters
H = 10;  //  semiheight of plate
L = H/2;  //  semiwidth of plate
a = 1; // semilength of crack (center of crack is at x1=x2=0)

// Discretization parameters
lc1 = .5; // element size at the border
lc2 = .5; // element size at the crack tip

// Domain construction
Point(11) = {0.0,0.0,0.0,lc1};
Point(12) = {L,0.0,0.0,lc1};
Point(21) = {L,2*H,0.0,lc1};
Point(22) = {0,2*H,0.0,lc1};
Point(3) = {L,H,0.0,lc1};
Point(4) = {0.0,H,0.0,lc1};

Line(11) = {11,12};
Line(12) = {12,3};
Line(13) = {3,4};
Line(14) = {4,11};
Line(21) = {21,22};
Line(23) = {4,3};
Line(22) = {22,4};
Line(24) = {3,21};

Line Loop(111) = {11,12,13,14};
Line Loop(211) = {21,22,23,24};
Plane Surface(1) = {111};
Plane Surface(2) = {211};

Physical Line(11) = {11};
Physical Line(12) = {12};
Physical Line(13) = {13};
Physical Line(23) = {23};
Physical Line(14) = {14};
Physical Line(21) = {21};
Physical Line(22) = {22};
Physical Line(24) = {24};

Physical Surface(5) = {1};
Physical Surface(6) = {2};

Geometry.Tolerance = 1e-5; // adjust value here for correct merge result
Coherence Mesh;Coherence;
Coherence;
