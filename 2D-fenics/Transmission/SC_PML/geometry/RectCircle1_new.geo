// Inputs

gridsize=1/13;
a=10;
b=10;
r=0.45;
L_pml=8;
xc=2;
Ty = 1;
Tx=1;
yc=5;


// Mesh definition
// Outer boundary
Point(1)={a+L_pml,-L_pml,0,gridsize};
Point(2)={a+L_pml,b+L_pml,0,gridsize};
Point(3)={0,b+L_pml,0,gridsize};
Point(4)={0,-L_pml,0,gridsize};
Point(7)={0,b,0,gridsize};
Point(8)={0,0,0,gridsize};

Point(992) = {a,0, 0.0, gridsize} ;
Point(993) = {a, b, 0.0, gridsize} ;


Point(9) = {xc, yc, 0.0, gridsize} ;
Point(10) = {xc, yc-r, 0.0, gridsize} ;
Point(11) = {xc-r, yc, 0.0, gridsize} ;
Point(12) = {xc, yc+r, 0.0, gridsize} ;
Point(13) = {xc+r, yc, 0.0, gridsize} ;
Circle(14) = {10, 9, 11};
Circle(15) = {11, 9, 12};
Circle(16) = {12, 9, 13};
Circle(17) = {13, 9, 10};

Translate {1*Tx, 0, 0} {
  Duplicata { Line{15, 16, 17, 14}; }
}
Translate {2*Tx, 0, 0} {
  Duplicata { Line{15, 16, 17, 14}; }
}
Translate {3*Tx, 0, 0} {
  Duplicata { Line{15, 16, 17, 14}; }
}
Translate {4*Tx, 0, 0} {
  Duplicata { Line{15, 16, 17, 14}; }
}


Line(1) = {7, 993};
Line(2) = {993, 992};
Line(3) = {992, 8};
Line(4) = {8, 7};
Line(5) = {7, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 4};
Line(9) = {4, 8};
Line Loop(34) = {1, 2, 3, -9, -8, -7, -6, -5};
Plane Surface(35) = {34};
Line Loop(36) = {1, 2, 3, 4};
Line Loop(37) = {17, 14, 15, 16};
Line Loop(38) = {21, 18, 19, 20};
Line Loop(39) = {25, 22, 23, 24};
Line Loop(40) = {29, 26, 27, 28};
Line Loop(41) = {33, 30, 31, 32};
Plane Surface(42) = {36, 37, 38, 39, 40, 41};

