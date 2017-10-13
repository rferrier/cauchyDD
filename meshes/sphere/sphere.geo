// see https://sites.google.com/site/auxcapucins/maillage-3d-en-gmsh---maillage-d-une-sphere
// Geometric parameters
R1 = 5;
R2 = 15;
pi = 3.1415926535;
// Discretization parameters

lc = 2;     // element size at the border

Point(1) = {0,0,0,lc};
Point(2) = {R1,0,0,lc};
Point(3) = {0,R1,0,lc};
Point(4) = {0,0,R1,lc};
Point(5) = {-R1,0,0,lc};
Point(6) = {0,-R1,0,lc};
Point(7) = {0,0,-R1,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,2};
Circle(5) = {2,1,7};
Circle(6) = {7,1,5};
Circle(7) = {5,1,4};
Circle(8) = {4,1,2};
Circle(9) = {6,1,7};
Circle(10) = {7,1,3};
Circle(11) = {3,1,4};
Circle(12) = {4,1,6};

Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11}; 
Line Loop(3) = {3,-12,-7}; 
Line Loop(4) = {4,-8,12}; 
Line Loop(5) = {5,10,-1}; 
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9}; 
Line Loop(8) = {-4,9,-5}; 

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

Surface Loop (1) = {1,2,3,4,5,6,7,8};
//Volume (1) = {1};
///////////////////////////////////////////

Point(102) = {R2,0,0,lc};
Point(103) = {0,R2,0,lc};
Point(104) = {0,0,R2,lc};
Point(105) = {-R2,0,0,lc};
Point(106) = {0,-R2,0,lc};
Point(107) = {0,0,-R2,lc};

Circle(101) = {102,1,103};
Circle(102) = {103,1,105};
Circle(103) = {105,1,106};
Circle(104) = {106,1,102};
Circle(105) = {102,1,107};
Circle(106) = {107,1,105};
Circle(107) = {105,1,104};
Circle(108) = {104,1,102};
Circle(109) = {106,1,107};
Circle(110) = {107,1,103};
Circle(111) = {103,1,104};
Circle(112) = {104,1,106};

Line Loop(101) = {101,111,108};
Line Loop(102) = {102,107,-111}; 
Line Loop(103) = {103,-112,-107}; 
Line Loop(104) = {104,-108,112}; 
Line Loop(105) = {105,110,-101}; 
Line Loop(106) = {-102,-110,106};
Line Loop(107) = {-103,-106,-109}; 
Line Loop(108) = {-104,109,-105}; 

Ruled Surface(101) = {101};
Ruled Surface(102) = {102};
Ruled Surface(103) = {103};
Ruled Surface(104) = {104};
Ruled Surface(105) = {105};
Ruled Surface(106) = {106};
Ruled Surface(107) = {107};
Ruled Surface(108) = {108};

Surface Loop (101) = {101,102,103,104,105,106,107,108};
Volume (1) = {101,1};

Physical Surface(1) = {101,102,103,104,105,106,107,108};
Physical Surface(2) = {1,2,3,4,5,6,7,8};
//Physical Surface(3) = {1}; // Loading surfaces
//Physical Surface(4) = {7}; //
Physical Volume(1) = {1};
