// Inputs

gridsize=1/30;
a=1;
theta = Pi/3;
R = 0.4;
radius_int = 0.3;
phi = Pi/6;
alpha = Asin((Sin(phi)*R)/radius_int);

xc = 0;
xcc = 0.75;
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
Point(5)={+xcc,ycc,0,gridsize};
Point(6)={R+xcc,ycc,0,gridsize};
Point(7)={-R+xcc,ycc,0,gridsize};
Point(8)={+xcc,-R+ycc,0,gridsize};
Point(9)={+xcc,+R+ycc,0,gridsize};

// inner circle resonator
Point(10)={xcc+radius_int,ycc,0,gridsize};
Point(11)={xcc-radius_int,ycc,0,gridsize};
Point(12)={xcc-radius_int,ycc,0,gridsize};
Point(13)={xcc+radius_int,ycc,0,gridsize};

Rotate {{0, 0, 1}, {xcc, ycc, 0}, phi} {
  Duplicata { Point{7}; }
}
Rotate {{0, 0, 1}, {xcc, ycc, 0}, -phi} {
  Duplicata { Point{7}; }
}


Point(99)={xcc-radius_int*Cos(alpha),ycc+radius_int*Sin(alpha),0,gridsize};
Point(100)={xcc-radius_int*Cos(alpha),ycc-radius_int*Sin(alpha),0,gridsize};



Circle(5) = {99, 5, 10};
Circle(6) = {10, 5, 100};
Line(7) = {100, 14};
Line(8) = {99, 15};
Circle(9) = {15, 5, 6};
Circle(10) = {6, 5, 14};

Periodic Line {1}={3};
Periodic Line {2}={4};

Delete {
  Point{11,7,9,8};
}

Line Loop(11) = {3, 4, 1, 2};Line Loop(12) = {9, 10, -7, -6, -5, 8};
Plane Surface(13) = {11, 12};
