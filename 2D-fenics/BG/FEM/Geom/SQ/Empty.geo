// Inputs

gridsize=1/10;
L = 1;
// Mesh definition
// Outer boundary
Point(1)={0.5*L,-0.5*L,0,gridsize};
Point(2)={0.5*L,0.5*L,0,gridsize};
Point(3)={-0.5*L,0.5*L,0,gridsize};
Point(4)={-0.5*L,-0.5*L,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(15) = {4, 1, 2, 3};
Periodic Line {1}={3};
Periodic Line {2}={4};
Plane Surface(6) = {15};
