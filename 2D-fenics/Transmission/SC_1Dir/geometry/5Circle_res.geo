// Inputs
gridsize=1/30;

//radius_int = 0.3;

radius = 0.45;

DefineConstant[ radius_int = {0.3, Name "R_in", Min 0.12, Max radius-0.05, Step 0.01} ];
a = 1;
xc=2;
yc=0.5;
phi = Pi/10;
alpha = Asin((Sin(phi)*radius)/radius_int);

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

// Inner circle
Point(5)={xc,yc,0,gridsize};
Point(6)={radius+xc,yc,0,gridsize};
Point(7)={xc-radius,yc,0,gridsize};
Point(8)={xc,-radius+yc,0,gridsize};
Point(9)={xc,radius+yc,0,gridsize};


// inner circle resonator
Point(10)={radius_int+xc,yc,0,gridsize};
Point(11)={-radius_int+xc,yc,0,gridsize};
Point(12)={xc,-radius_int+yc,0,gridsize};
Point(13)={xc,radius_int+yc,0,gridsize};


Point(99)={-radius_int*Cos(alpha)+xc,radius_int*Sin(alpha)+yc,0,gridsize};
Point(100)={-radius_int*Cos(alpha)+xc,yc-radius_int*Sin(alpha),0,gridsize};
Point(17)={-radius*Cos(phi)+xc,radius*Sin(phi)+yc,0,gridsize};
Point(16)={-radius*Cos(phi)+xc,-radius*Sin(phi)+yc,0,gridsize};


Circle(5) = {6, 5, 17};
Circle(6) = {6, 5, 16};
Circle(7) = {99, 5, 10};
Circle(8) = {10, 5, 100};
Line(9) = {17, 99};
Line(10) = {16, 100};

Translate {a, 0, 0} {
  Duplicata { Line{7, 8, 10, 6, 5, 9}; }
}
Translate {2*a, 0, 0} {
  Duplicata { Line{7, 8, 10, 6, 5, 9}; }
}
Translate {3*a, 0, 0} {
  Duplicata { Line{7, 8, 10, 6, 5, 9}; }
}
Translate {4*a, 0, 0} {
  Duplicata { Line{7, 8, 10, 6, 5, 9}; }
}

Periodic Line {1}={3};
Periodic Line {2}={4};
Line Loop(35) = {2, 3, 4, 1};
Line Loop(36) = {7, 8, -10, -6, 5, 9};
Line Loop(37) = {11, 12, -13, -14, 15, 16};
Line Loop(38) = {17, 18, -19, -20, 21, 22};
Line Loop(39) = {23, 24, -25, -26, 27, 28};
Line Loop(40) = {29, 30, -31, -32, 33, 34};
Plane Surface(41) = {35, 36, 37, 38, 39, 40};
