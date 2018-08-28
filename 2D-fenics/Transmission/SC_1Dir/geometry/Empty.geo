// Inputs
gridsize=1/30;
radius=0.45;
a = 1;
xc=2;
yc=0.5;


// Mesh definition
// Outer boundary
Point(1)={8*a,0,0,gridsize};
Point(2)={8*a,a,0,gridsize};
Point(3)={0,a,0,gridsize};
Point(4)={0,0,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};




Line Loop(28) = {4, 1, 2, 3};

Periodic Line {1}={3};
Periodic Line {2}={4};
Plane Surface(30) = {28};

