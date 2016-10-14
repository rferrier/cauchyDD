// Enable cracks (or not...)
Geometry.AutoCoherence = 0;

// Geometric parameters
pi    = 3.1415926535;  // What I love to make learn a number useful to the wise
R     = 2;
a1    = 2.5;   // Radius of the crack
b1    = 1.5;

p     = b1*b1/a1;  // Ellipse parameters
e     = Sqrt(a1*a1-b1*b1)/a1;

a     = 4;      // Center of the crack
b     = 3;      //
H     = 10;  //  height of plate
L     = 7;  //  width of plate

// Discretization parameters
lc1 = 1; // element size at the border
lc2 = 1; // element size at the crack tip

// Domain construction
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L,0.0,0.0,lc1};
Point(3) = {L,L,0.0,lc2};
Point(4) = {0.0,L,0.0,lc2};
Point(5) = {0.0,0.0,H,lc1};
Point(6) = {L,0.0,H,lc1};
Point(7) = {L,L,H,lc2};
Point(8) = {0.0,L,H,lc2};

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

Line Loop(101) = {1,2,3,4};
//Plane Surface(1) = {101,107}; (see under it)
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
Point(9) = {a1+a,b,0,lc1};
Point(10) = {b1*Cos(pi/2)+a,b1*Sin(pi/2)+b,0,lc1};
Point(11) = {a1*Cos(pi)+a,a1*Sin(pi)+b,0,lc1};
Point(12) = {b1*Cos(3*pi/2)+a,b1*Sin(3*pi/2)+b,0,lc1};
Point(13) = {a,b,0,lc1};

Circle(13) = {9,13,9,10};
Circle(14) = {10,13,9,11};
Circle(15) = {11,13,9,12};
Circle(16) = {12,13,9,9};

//Point(9) = {R+a,b,0,lc1};
//Point(10) = {R*Cos(pi/2)+a,R*Sin(pi/2)+b,0,lc1};
//Point(11) = {R*Cos(pi)+a,R*Sin(pi)+b,0,lc1};
//Point(12) = {R*Cos(3*pi/2)+a,R*Sin(3*pi/2)+b,0,lc1};
//Point(13) = {a,b,0,lc1};

//Circle(13) = {9,13,10};
//Circle(14) = {10,13,11};
//Circle(15) = {11,13,12};
//Circle(16) = {12,13,9};

Line Loop(107) = {13,14,15,16};
Plane Surface(1) = {101,107};
Plane Surface(7) = {107};

// Physical stuff
Surface Loop(1001) = {1,2,3,4,5,6,7};

Volume(9) = {1001};


Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
//Physical Surface(70) = {-7};

Physical Volume(9) = {9};

