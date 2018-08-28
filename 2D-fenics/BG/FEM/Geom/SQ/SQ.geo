// Inputs

gridsize=1/30;
d=0.7;

// Mesh definition
// Outer boundary
Point(1)={0.5,-0.5,0,gridsize};
Point(2)={0.5,0.5,0,gridsize};
Point(3)={-0.5,0.5,0,gridsize};
Point(4)={-0.5,-0.5,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Inner circle
Point(6)={d/2,-d/2,0,gridsize};
Point(7)={d/2,d/2,0,gridsize};
Point(8)={-d/2,d/2,0,gridsize};
Point(9)={-d/2,-d/2,0,gridsize};

Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,9};
Line(13) = {9,6};

Line Loop(14) = {13, 10, 11, 12};

Line Loop(15) = {4, 1, 2, 3};
Periodic Line {1}={3};
Periodic Line {2}={4};
Plane Surface(6) = {15,14};
