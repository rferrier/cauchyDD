
// Geometric parameters
H = 10;  //  height of plate
L = 5;  //  width of plate

// Discretization parameters
lc1 = .5; // element size at the border
lc2 = .5; // element size at the crack tip

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,H,0.0,lc2};
Point(4) = {0.0,H,0.0,lc2};

// Subdomain limits
Point(5)  = {0.0,H/5,0.0,lc1};
Point(6)  = {L,H/5,0.0,lc1};
Point(7)  = {0.0,2*H/5,0.0,lc1};
Point(8)  = {L,2*H/5,0.0,lc1};
Point(9)  = {0.0,3*H/5,0.0,lc1};
Point(10) = {L,3*H/5,0.0,lc1};
Point(11) = {0.0,4*H/5,0.0,lc1};
Point(12) = {L,4*H/5,0.0,lc1};

Line(13) = {5,6};
Line(14) = {7,8};
Line(15) = {9,10};
Line(16) = {11,12};

Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,8};
Line(4) = {8,10};
Line(5) = {10,12};
Line(6) = {12,3};
Line(7) = {3,4};
Line(8) = {4,11};
Line(9) = {11,9};
Line(10) = {9,7};
Line(11) = {7,5};
Line(12) = {5,1};

Line Loop(11) = {1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(1) = {11};

// Add elements
Line{13} In Surface{1};
Line{14} In Surface{1};
Line{15} In Surface{1};
Line{16} In Surface{1};

//
Physical Line(1) = {1};
Physical Line(2) = {2,3,4,5,6};
Physical Line(3) = {7};
Physical Line(4) = {8,9,10,11,12}; 
Physical Surface(5) = {1};
