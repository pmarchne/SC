// Inputs

gridsize=1/50;
a=1;
theta = Pi/3;
R = 0.4;

xc = 0;
xcc = -0.75;
yc = 0;
ycc = Sin(theta)/2;

// Mesh definition
// Outer boundary
Point(1)={a-xc,-yc,0,gridsize};
Point(2)={a*(Cos(theta)+1)-xc,a*Sin(theta)-yc,0,gridsize};
Point(3)={a*Cos(theta)-xc,a*Sin(theta)-yc,0,gridsize};
Point(4)={-xc,-yc,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Inner circle
Point(5)={-xcc,ycc,0,gridsize};
Point(6)={R-xcc,ycc,0,gridsize};
Point(7)={-R-xcc,ycc,0,gridsize};
Point(8)={-xcc,-R+ycc,0,gridsize};
Point(9)={-xcc,+R+ycc,0,gridsize};

Circle(5) = {9, 5, 6};
Circle(6) = {6, 5, 8};
Circle(7) = {8, 5, 7};
Circle(8) = {7, 5, 9};

Line Loop(9) = {5, 6, 7, 8};

Line Loop(10) = {4, 1, 2, 3};
Periodic Line {1}={3};
Periodic Line {2}={4};
Plane Surface(6) = {10,9};
