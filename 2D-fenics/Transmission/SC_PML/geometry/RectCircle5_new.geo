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


Translate {0, Ty, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, -Ty, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, 2*Ty, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, -2*Ty, 0} {
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

Line Loop(114) = {1, 2, 3, -9, -8, -7, -6, -5};
Plane Surface(115) = {114};
Line Loop(116) = {4, 1, 2, 3};
Line Loop(117) = {44, 45, 42, 43};
Line Loop(118) = {30, 31, 32, 33};
Line Loop(119) = {34, 38, 36, 37};
Line Loop(120) = {35, 39, 40, 41};
Line Loop(121) = {29, 26, 27, 28};
Line Loop(122) = {46, 47, 48, 49};
Line Loop(123) = {50, 51, 52, 53};
Line Loop(124) = {54, 58, 56, 57};
Line Loop(125) = {55, 59, 60, 61};
Line Loop(126) = {15, 16, 17, 14};
Line Loop(127) = {18, 19, 20, 21};
Line Loop(128) = {22, 23, 24, 25};
Line Loop(129) = {100, 101, 95, 99};
Line Loop(130) = {87, 88, 89, 86};
Line Loop(131) = {90, 91, 92, 93};
Line Loop(132) = {94, 98, 96, 97};
Line Loop(133) = {102, 103, 104, 105};
Line Loop(134) = {106, 107, 108, 109};
Line Loop(135) = {110, 111, 112, 113};
Line Loop(136) = {72, 73, 70, 71};
Line Loop(137) = {62, 63, 64, 65};
Line Loop(138) = {66, 67, 68, 69};
Line Loop(139) = {74, 78, 76, 77};
Line Loop(140) = {75, 79, 80, 81};
Line Loop(141) = {82, 83, 84, 85};
Plane Surface(142) = {116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141};
