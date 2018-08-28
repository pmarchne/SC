// Inputs
radius=0.1;
radius_int = 0.05;
a = 0.3;
gridsize=a/40;
phi = Asin(0.5/2);
alpha = Asin((Sin(phi)*radius)/radius_int);

// Mesh definition
// Outer boundary
Point(1)={a/2,-a/2,0,gridsize};
Point(2)={a/2,a/2,0,gridsize};
Point(3)={-a/2,a/2,0,gridsize};
Point(4)={-a/2,-a/2,0,gridsize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Inner circle
Point(5)={0.0,0.0,0,gridsize};
Point(6)={radius,0.0,0,gridsize};
Point(7)={-radius,0.0,0,gridsize};
Point(8)={0.0,-radius,0,gridsize};
Point(9)={0.0,radius,0,gridsize};


// inner circle resonator
Point(10)={radius_int,0.0,0,gridsize};
Point(11)={-radius_int,0.0,0,gridsize};
Point(12)={0.0,-radius_int,0,gridsize};
Point(13)={0.0,radius_int,0,gridsize};

Rotate {{0, 0, 1}, {0, 0, 0}, phi} {
  Duplicata { Point{11}; }
}
Rotate {{0, 0, 1}, {0, 0, 0}, -phi} {
  Duplicata { Point{11}; }
}



Rotate {{0, 0, 1}, {0, 0, 0}, phi} {
  Duplicata { Point{7}; }
}
Rotate {{0, 0, 1}, {0, 0, 0}, -phi} {
  Duplicata { Point{7}; }
}


Point(99)={-radius_int*Cos(alpha),radius_int*Sin(alpha),0,gridsize};
Point(100)={-radius_int*Cos(alpha),-radius_int*Sin(alpha),0,gridsize};

Circle(5) = {6, 5, 17};
Circle(6) = {6, 5, 16};
Circle(7) = {99, 5, 10};
Circle(8) = {10, 5, 100};
Line(9) = {17, 99};
Line(10) = {16, 100};

Periodic Line {1}={3};
Periodic Line {2}={4};Line Loop(11) = {3, 4, 1, 2};
Line Loop(12) = {5, 9, 7, 8, -10, -6};
Plane Surface(13) = {11, 12};
Delete {
  Point{15, 14, 11};
}
Rotate {{0, 0, 1}, {0, 0, 0}, 0} {
  Line{5, 9, 7, 8, 10, 6};
}
