// Enable cracks (or not...)
//Geometry.AutoCoherence = 0;

// Geometric parameters
pi    = 3.1415926535;  // What I love to make learn a number useful to wises
a1    = 2;   // Elliptic crack
b1    = 1;
a2    = 1.5;   // Elliptic crack
b2    = 2.5;
a     = 4;      // Centers of the cracks
b     = 2;      //
c     = 4;      //
d     = 7;      //
da    = pi/15; // angle of the cracks
H     = 2;  //  height of plate 2-6
L     = 7;  //  width of plate
E     = 10; //E     = 7;
eps   = .4;  //refining width

// Discretization parameters
lc1 = .5; // element size at the border
lc2 = .5; // element size at the crack tip
lc3 = .5; // refining size

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,E,0.0,lc2};
Point(4) = {0.0,E,0.0,lc2};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,E,H,lc2};
Point(8) = {0.0,E,H,lc2};

// Refining points
Point(101) = {eps,eps,eps,lc2};
Point(102) = {L-eps,eps,eps,lc2};
Point(103) = {L-eps,E-eps,eps,lc2};
Point(104) = {eps,E-eps,eps,lc2};
Point(105) = {eps,eps,H-eps,lc2};
Point(106) = {L-eps,eps,H-eps,lc2};
Point(107) = {L-eps,E-eps,H-eps,lc2};
Point(108) = {eps,E-eps,H-eps,lc2};

// Cut points
Point(21)  = {0,0,H/2-a*Tan(da),lc1};
Point(22)  = {L,0,H/2+(L-a)*Tan(da),lc1};
Point(23)  = {L,E,H/2+(L-a)*Tan(da),lc1};
Point(24)  = {0,E,H/2-a*Tan(da),lc1};
Point(121) = {a-(H/2-eps)/Tan(da),eps,eps,lc1};
Point(122) = {L-eps,eps,H/2+(L-a-eps)*Tan(da),lc1};
Point(123) = {L-eps,E-eps,H/2+(L-a-eps)*Tan(da),lc1};
Point(124) = {a-(H/2-eps)/Tan(da),E-eps,eps,lc1};

// Aretes
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {6,5};
Line(6) = {7,6};
Line(7) = {8,7};
Line(8) = {5,8};

Line(9) = {1,21};
Line(2009) = {21,5};
Line(10) = {2,22};
Line(2010) = {22,6};
Line(11) = {3,23};
Line(2011) = {23,7};
Line(12) = {4,24};
Line(2012) = {24,8};

// refining aretes
Line(101) = {101,121}; //121
Line(2101) = {121,102};
Line(102) = {102,103};
Line(103) = {103,124}; //124
Line(2103) = {124,104};
Line(104) = {104,101};

Line(105) = {106,105};
Line(106) = {107,106};
Line(107) = {108,107};
Line(108) = {105,108};

Line(109) = {101,105};
Line(110) = {102,122};
Line(2110) = {122,106};
Line(111) = {103,123};
Line(2111) = {123,107};
Line(112) = {104,108};

//Cut lines
Line(201) = {21,22};
Line(202) = {22,23};
Line(203) = {23,24};
Line(204) = {24,21};
Line(211) = {121,122};
Line(212) = {122,123};
Line(213) = {123,124};
Line(214) = {124,121};

// External faces
Line Loop(101) = {1,2,3,4};
Plane Surface(1) = {101};
Line Loop(102) = {5,6,7,8};
Plane Surface(2) = {102};
Line Loop(103) = {1,10,2010,5,-2009,-9};
Plane Surface(3) = {103};
Line Loop(104) = {2,11,2011,6,-2010,-10};
Plane Surface(4) = {104};
Line Loop(105) = {3,12,2012,7,-2011,-11};
Plane Surface(5) = {105};
Line Loop(106) = {4,9,2009,8,-2012,-12};
Plane Surface(6) = {106};

// Refining faces
Line Loop(1101) = {101,2101,102,103,2103,104};
Plane Surface(101) = {1101};
Line Loop(1102) = {105,106,107,108};
Plane Surface(102) = {1102};
Line Loop(1103) = {101,2101,110,2110,105,-109};
Plane Surface(103) = {1103};
Line Loop(1104) = {102,111,2111,106,-2110,-110};
Plane Surface(104) = {1104};
Line Loop(1105) = {103,2103,112,107,-2111,-111};
Plane Surface(105) = {1105};
Line Loop(1106) = {104,109,108,-112};
Plane Surface(106) = {1106};

// Crack construction
Point(9) = {a,b,H/2,lc1};
Point(10) = {a,b+b1,H/2,lc1};
Point(11) = {a-a1*Cos(da),b,-a1*Sin(da)+H/2,lc1};
Point(12) = {a,b-b1,H/2,lc1};
Point(13) = {a+a1*Cos(da),b,a1*Sin(da)+H/2,lc1};

Ellipse(13) = {10,9,11,11};  //{start,center,major,end}
Ellipse(14) = {11,9,11,12};
Ellipse(15) = {12,9,11,13};
Ellipse(16) = {13,9,11,10};

Point(14) = {c,d,H/2,lc1};
Point(15) = {c,d+b2,H/2,lc1};
Point(16) = {c-a2*Cos(da),d,-a2*Sin(da)+H/2,lc1};
Point(17) = {c,d-b2,H/2,lc1};
Point(18) = {c+a2*Cos(da),d,a2*Sin(da)+H/2,lc1};

Ellipse(17) = {15,14,16,16};  //{start,center,major,end}
Ellipse(18) = {16,14,16,17};
Ellipse(19) = {17,14,16,18};
Ellipse(20) = {18,14,16,15};

Line Loop(107) = {13,14,15,16};
Plane Surface(7) = {107};
Line Loop(108) = {17,18,19,20};
Plane Surface(8) = {108};

//Cut faces
Line Loop(201) = {201,202,203,204};//,-211,-212,-213,-214};
Line Loop(202) = {211,212,213,214};//,-13,-14,-15,-16,-17,-18,-19,-20};
Plane Surface(201) = {201,202};
Plane Surface(202) = {202,107,108};

// Physical stuff
Surface Loop(1001) = {1,2,3,4,5,6};
Volume(9) = {1001};

/*Point{21} In Surface{3};
Point{21} In Surface{6};
Point{22} In Surface{4};
Point{22} In Surface{3};
Point{23} In Surface{4};
Point{23} In Surface{5};
Point{24} In Surface{5};
Point{24} In Surface{6};

Point{121} In Surface{103};
Point{121} In Surface{101};
Point{122} In Surface{104};
Point{122} In Surface{103};
Point{123} In Surface{104};
Point{123} In Surface{105};
Point{124} In Surface{105};
Point{124} In Surface{101};*/

Line{214} In Surface{101};
Line{213} In Surface{105};
Line{212} In Surface{104};
Line{211} In Surface{103};

Line{204} In Surface{6};
Line{203} In Surface{5};
Line{202} In Surface{4};
Line{201} In Surface{3};

Surface{7} In Volume{9};
Surface{8} In Volume{9};
Surface{201} In Volume{9};
Surface{202} In Volume{9};
Surface{101} In Volume{9};
Surface{102} In Volume{9};
Surface{103} In Volume{9};
Surface{104} In Volume{9};
Surface{105} In Volume{9};
Surface{106} In Volume{9};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7,8};
Physical Surface(8) = {7,8,201,202};

Physical Volume(9) = {9};
