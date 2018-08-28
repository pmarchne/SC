// Inputs

gridsize=1/13;
a=10;
b=10;
radius=0.45;
L_pml=8;
radius_int = 0.4;
xc=2;
Ty = 1;
Tx=1;
yc=5;
phi = Pi/5;
alpha = Asin((Sin(phi)*radius)/radius_int);



// Mesh definition
// Outer boundary
Point(1)={a+L_pml,-L_pml,0,gridsize};
Point(2)={a+L_pml,b+L_pml,0,gridsize};
Point(3)={0,b+L_pml,0,gridsize};
Point(4)={0,-L_pml,0,gridsize};
Point(997)={0,b,0,gridsize};
Point(998)={0,0,0,gridsize};

Point(992) = {a,0, 0.0, gridsize} ;
Point(993) = {a, b, 0.0, gridsize} ;


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


Circle(55) = {6, 5, 17};
Circle(66) = {6, 5, 16};
Circle(77) = {99, 5, 10};
Circle(88) = {10, 5, 100};

Line(1) = {997, 993};
Line(2) = {993, 992};
Line(3) = {992, 998};
Line(4) = {998, 997};
Line(5) = {997, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 4};
Line(9) = {4, 998};
Line(11) = {17, 99};
Line(12) = {16, 100};

Line Loop(194) = {5, 6, 7, 8, 9, -3, -2, -1};
Plane Surface(195) = {194};
Delete {
  Point{7, 11};
}

Translate {Tx, 0, 0} {
  Duplicata { Line{77, 55, 11, 12, 66, 88}; }
}
Translate {2*Tx, 0, 0} {
  Duplicata { Line{77, 55, 11, 12, 66, 88}; }
}
Translate {3*Tx, 0, 0} {
  Duplicata { Line{77, 55, 11, 12, 66, 88}; }
}
Translate {4*Tx, 0, 0} {
  Duplicata { Line{77, 55, 11, 12, 66, 88}; }
}
Line Loop(220) = {1, 2, 3, 4};
Line Loop(221) = {77, 88, -12, -66, 55, 11};
Line Loop(222) = {202, 207, -205, -206, 203, 204};
Line Loop(223) = {201, -199, -200, 197, 198, 196};
Line Loop(224) = {211, -213, -208, -210, -209, 212};
Line Loop(225) = {219, -217, -218, 215, 216, 214};
Plane Surface(226) = {220, 221, 222, 223, 224, 225};

