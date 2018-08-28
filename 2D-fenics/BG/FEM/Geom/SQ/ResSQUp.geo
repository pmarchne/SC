// Inputs

gridsize=1/15;
ratio=1;
d=0.7*ratio;
g=0.1*ratio;
a=0.5*ratio;

// Mesh definition
// Outer boundary
Point(1)={a,-a,0,gridsize};
Point(2)={a,a,0,gridsize};
Point(3)={-a,a,0,gridsize};
Point(4)={-a,-a,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Inner circle
Point(6)={d/2,-d/2,0,gridsize};
Point(7)={d/2,d/2,0,gridsize};
Point(8)={-d/2,d/2,0,gridsize};
Point(9)={-d/2,-d/2,0,gridsize};

Point(60)={d/2-g,-d/2+g,0,gridsize};
Point(70)={d/2-g,d/2,0,gridsize};
Point(80)={-d/2+g,d/2,0,gridsize};
Point(90)={-d/2+g,-d/2+g,0,gridsize};


Line(10) = {6,7};
Line(11) = {7,70};
Line(12) = {70,60};
Line(13) = {60,90};
Line(14) = {90,80};
Line(15) = {80,8};
Line(16) = {8,9};
Line(17) = {9,6};

Line Loop(24) = {17,10,11,12,13,14,15,16};
Line Loop(25) = {4, 1, 2, 3};

Periodic Line {1}={3};
Periodic Line {2}={4};
Plane Surface(6) = {25,24};