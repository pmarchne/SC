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


Translate {0, 1, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, -1, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
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

Line Loop(74) = {1, 2, 3, -9, -8, -7, -6, -5};
Plane Surface(75) = {74};
Line Loop(76) = {4, 1, 2, 3};
Line Loop(77) = {32, 33, 30, 31};
Line Loop(78) = {24, 25, 22, 23};
Line Loop(79) = {26, 27, 28, 29};
Line Loop(80) = {34, 38, 36, 37};
Line Loop(81) = {35, 39, 40, 41};
Line Loop(82) = {15, 16, 17, 14};
Line Loop(83) = {18, 19, 20, 21};
Line Loop(84) = {66, 67, 68, 69};
Line Loop(85) = {58, 56, 57, 54};
Line Loop(86) = {59, 60, 61, 55};
Line Loop(87) = {62, 63, 64, 65};
Line Loop(88) = {70, 71, 72, 73};
Line Loop(89) = {49, 46, 47, 48};
Line Loop(90) = {42, 43, 44, 45};
Line Loop(91) = {50, 51, 52, 53};
Plane Surface(92) = {76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91};
