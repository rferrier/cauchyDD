
// Geometric parameters
H   = 1;     //  height of plate
L   = 1;     //  width of plate
E   = 1;
a   = .3;     // Center of the crack
b   = .2;     //
c   = .5;
d   = .7;
a1  = .15;   // Radius of the crack
b1  = .1;
a2  = .15;   // Radius of the crack
b2  = .2;
// Discretization parameters

lc1 = .2;     // element size at the border

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

// Ellipse 1
Point(9) = {a,b,H,lc1};
Point(10) = {a,b+b1,H,lc1};
Point(11) = {a-a1,b,H,lc1};
Point(12) = {a,b-b1,H,lc1};
Point(13) = {a+a1,b,H,lc1};

Ellipse(13) = {10,9,11,11};  //{start,center,major,end}
Ellipse(14) = {11,9,11,12};
Ellipse(15) = {12,9,11,13};
Ellipse(16) = {13,9,11,10};

// Ellipse 2
Point(14) = {c,d,H,lc1};
Point(15) = {c,d+b2,H,lc1};
Point(16) = {c-a2,d,H,lc1};
Point(17) = {c,d-b2,H,lc1};
Point(18) = {c+a2,d,H,lc1};

Ellipse(17) = {15,14,16,16};  //{start,center,major,end}
Ellipse(18) = {16,14,16,17};
Ellipse(19) = {17,14,16,18};
Ellipse(20) = {18,14,16,15};

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
Line Loop(108) = {17,18,19,20};
Plane Surface(2) = {102,107,108};
Plane Surface(7) = {107};
Plane Surface(8) = {108};

Line{13,14,15,16} In Surface{2};
Line{17,18,19,20} In Surface{2};

Surface Loop(1001) = {1,2,3,4,5,6,7,8};
Volume(1) = {1001};

Physical Surface(1) = {1};
Physical Surface(2) = {2,7,8};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7,8};

Physical Volume(1) = {1};
