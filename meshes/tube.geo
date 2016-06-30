
// Geometric parameters
pi    = 3.1415926535;
R     = 10;    // interior radius
Re    = 15;    // exterior radius
alpha = pi/10;  // encastrement angle
theta = pi/6;  // fluid level

// Discretization parameters
lc1 = 1; // element size in the interior
lc2 = 1; // element size at exterior

// Domain construction
Point(1) = {0.0,-Re,0.0,lc2};
Point(2) = {Re*Cos(pi/2-alpha),Re*Sin(pi/2-alpha),0.0,lc2};
Point(3) = {Re*Cos(pi/2+alpha),Re*Sin(pi/2+alpha),0.0,lc2};
Point(4) = {0.0,-R,0.0,lc1};
Point(5) = {R*Cos(theta),R*Sin(theta),0.0,lc1};
Point(6) = {R*Cos(pi-theta),R*Sin(pi-theta),0.0,lc1};
Point(7) = {0.0,0.0,0.0,lc1};

Circle(1) = {1,7,2};
Circle(2) = {2,7,3};
Circle(3) = {3,7,1};
Circle(4) = {4,7,5};
Circle(5) = {5,7,6};
Circle(6) = {6,7,4};

Line Loop(11) = {1,2,3};
Line Loop(12) = {4,5,6};
Plane Surface(1) = {11,12};

Physical Line(1)    = {1,3};
Physical Line(2)    = {2};
Physical Line(3)    = {4,5,6};
Physical Line(8)    = {4,6}; 
Physical Line(9)    = {5};
Physical Surface(4) = {1};
