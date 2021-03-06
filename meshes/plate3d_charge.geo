
// Geometric parameters
H   = 2;     //  height of plate
L   = 7;     //  width of plate
E   = 7;
a   = 4;     // Center of the crack
b   = 2;     //
a1  = 1.5;   // Radius of the crack
b1  = 1;
// Discretization parameters

lc1 = .75;     // element size at the border

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,E,0.0,lc1};
Point(4) = {0.0,E,0.0,lc1};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,E,H,lc1};
Point(8) = {0.0,E,H,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {6,5};
Line(6) = {7,6};
Line(7) = {8,7};
Line(8) = {5,8};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// Ellipse
Point(9) = {a,b,H,lc1};
Point(10) = {a,b+b1,H,lc1};
Point(11) = {a-a1,b,H,lc1};
Point(12) = {a,b-b1,H,lc1};
Point(13) = {a+a1,b,H,lc1};

Ellipse(13) = {10,9,11,11};  //{start,center,major,end}
Ellipse(14) = {11,9,11,12};
Ellipse(15) = {12,9,11,13};
Ellipse(16) = {13,9,11,10};

Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Line Loop(103) = {1,10,5,-9};
Plane Surface(3) = {103};
Line Loop(104) = {2,11,6,-10};
Plane Surface(4) = {104};
Line Loop(105) = {3,12,7,-11};
Plane Surface(5) = {105};
Line Loop(106) = {4,9,8,-12};
Plane Surface(6) = {106};
Line Loop(107) = {13,14,15,16};
Plane Surface(2) = {102,107};
Plane Surface(7) = {107};

Line{13,14,15,16} In Surface{2};

Surface Loop(1001) = {1,2,3,4,5,6,7};
Volume(1) = {1001};

Physical Surface(1) = {1};
Physical Surface(2) = {2,7};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};

Physical Volume(1) = {1};
