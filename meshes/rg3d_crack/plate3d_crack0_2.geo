// Enable cracks (or not...)
Geometry.AutoCoherence = 0;

// Geometric parameters
pi    = 3.1415926535;  // What I love to make learn a number useful to wises
a1    = .25;   // Elliptic crack
b1    = .15;
a     = .4;      // Centers of the cracks
b     = .5;      //
da    = pi/6; // angle of the cracks
db    = pi/10;
H     = 1;  //  height of plate 2-6
L     = 1;  //  width of plate
E     = 1; //E     = 7;

// Discretization parameters
lc1 = .065; // element size at the border
lc2 = .065; // element size at the crack tip
lc3 = .065; // refining size

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,E,0.0,lc2};
Point(4) = {0.0,E,0.0,lc2};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,E,H,lc2};
Point(8) = {0.0,E,H,lc2};

// Aretes
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

// External faces
Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Plane Surface(2) = {102};
Line Loop(103) = {1,10,5,-9};
Plane Surface(3) = {103};
Line Loop(104) = {2,11,6,-10};
Plane Surface(4) = {104};
Line Loop(105) = {3,12,7,-11};
Plane Surface(5) = {105};
Line Loop(106) = {4,9,8,-12};
Plane Surface(6) = {106};

// Crack construction
Point(9) = {a,b,H/2,lc1};
Point(10) = {a,b+b1*Cos(db),H/2+b1*Sin(db),lc1};
Point(11) = {a-a1*Cos(da),b,-a1*Sin(da)+H/2,lc1};
Point(12) = {a,b-b1*Cos(db),H/2-b1*Sin(db),lc1};
Point(13) = {a+a1*Cos(da),b,a1*Sin(da)+H/2,lc1};

Ellipse(13) = {10,9,11,11};  //{start,center,major,end}
Ellipse(14) = {11,9,11,12};
Ellipse(15) = {12,9,11,13};
Ellipse(16) = {13,9,11,10};

Line Loop(107) = {13,14,15,16};
Plane Surface(7) = {107};

// Physical stuff
Surface Loop(1001) = {1,2,3,4,5,6};
Volume(9) = {1001};

Surface{7} In Volume{9};
//Surface{8} In Volume{9};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};

Physical Volume(9) = {9};
