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
lc1 = .025; // fine element
lc2 = .0125; // element size at the crack tip
lc3 = .04; //  coarse size
lc4 = .05; //coarser size

// Domain construction
Point(1) = {0.0,0.0,0.0,lc3};
Point(2) = {L,0.0,0.0,lc3};
Point(3) = {L,E,0.0,lc3};
Point(4) = {0.0,E,0.0,lc3};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,E,H,lc1};
Point(8) = {0.0,E,H,lc1};

Point(19)  = {0.0,0.0,H/2,lc3};
Point(20) = {L,0.0,H/2,lc3};
Point(21) = {L,E,H/2,lc3};
Point(22) = {0.0,E,H/2,lc3};

Point(111) = {lc4,lc4,0.0,lc4};
Point(121) = {L-lc4,lc4,0.0,lc4};
Point(131) = {L-lc4,E-lc4,0.0,lc4};
Point(141) = {lc4,E-lc4,0.0,lc4};
Point(151) = {lc4,lc4,H,lc4};
Point(161) = {L-lc4,lc4,H,lc4};
Point(171) = {L-lc4,E-lc4,H,lc4};
Point(181) = {lc4,E-lc4,H,lc4};

// Aretes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {6,5};
Line(6) = {7,6};
Line(7) = {8,7};
Line(8) = {5,8};

Line(11) = {121,111};
Line(21) = {131,121};
Line(31) = {141,131};
Line(41) = {111,141};

Line(51) = {161,151};
Line(61) = {171,161};
Line(71) = {181,171};
Line(81) = {151,181};

Line(91)  = {1,19};
Line(101) = {2,20};
Line(111) = {3,21};
Line(121) = {4,22};
Line(92)  = {19,5};
Line(102) = {20,6};
Line(112) = {21,7};
Line(122) = {22,8};

Line(113) = {19,20};
Line(114) = {20,21};
Line(115) = {21,22};
Line(116) = {22,19};

// External faces
Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Plane Surface(2) = {102};
Line Loop(103) = {1,101,102,5,-92,-91};
Plane Surface(3) = {103};
Line Loop(104) = {2,111,112,6,-102,-101};
Plane Surface(4) = {104};
Line Loop(105) = {3,121,122,7,-112,-111};
Plane Surface(5) = {105};
Line Loop(106) = {4,91,92,8,-122,-121};
Plane Surface(6) = {106};

// Crack construction
Point(9) = {a,b,H/2,lc2};
Point(10) = {a,b+b1*Cos(db),H/2+b1*Sin(db),lc2};
Point(11) = {a-a1*Cos(da),b,-a1*Sin(da)+H/2,lc2};
Point(12) = {a,b-b1*Cos(db),H/2-b1*Sin(db),lc2};
Point(13) = {a+a1*Cos(da),b,a1*Sin(da)+H/2,lc2};

Ellipse(13) = {10,9,11,11};  //{start,center,major,end}
Ellipse(14) = {11,9,11,12};
Ellipse(15) = {12,9,11,13};
Ellipse(16) = {13,9,11,10};

Line Loop(107) = {13,14,15,16};
Plane Surface(7) = {107};

Line{113} In Surface{3};
Line{114} In Surface{4};
Line{115} In Surface{5};
Line{116} In Surface{6};

Line{11} In Surface{1};
Line{21} In Surface{1};
Line{31} In Surface{1};
Line{41} In Surface{1};

Line{51} In Surface{2};
Line{61} In Surface{2};
Line{71} In Surface{2};
Line{81} In Surface{2};

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
