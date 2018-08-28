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
Translate {0, 3*Ty, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, -3*Ty, 0} {
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
Line(9) = {4, 8};Line Loop(154) = {5, 6, 7, 8, 9, -3, -2, -1};
Plane Surface(155) = {154};
Line Loop(156) = {4, 1, 2, 3};
Line Loop(157) = {45, 42, 43, 44};
Line Loop(158) = {54, 58, 56, 57};
Line Loop(159) = {53, 50, 51, 52};
Line Loop(160) = {49, 46, 47, 48};
Line Loop(161) = {55, 59, 60, 61};
Line Loop(162) = {41, 35, 39, 40};
Line Loop(163) = {38, 36, 37, 34};
Line Loop(164) = {65, 62, 63, 64};
Line Loop(165) = {74, 78, 76, 77};
Line Loop(166) = {73, 70, 71, 72};
Line Loop(167) = {69, 66, 67, 68};
Line Loop(168) = {15, 16, 17, 14};
Line Loop(169) = {25, 22, 23, 24};
Line Loop(170) = {33, 30, 31, 32};
Line Loop(171) = {29, 26, 27, 28};
Line Loop(172) = {75, 79, 80, 81};
Line Loop(173) = {21, 18, 19, 20};
Line Loop(174) = {124, 125, 122, 123};
Line Loop(175) = {133, 130, 131, 132};
Line Loop(176) = {129, 126, 127, 128};
Line Loop(177) = {134, 138, 136, 137};
Line Loop(178) = {121, 115, 119, 120};
Line Loop(179) = {118, 116, 117, 114};
Line Loop(180) = {144, 145, 142, 143};
Line Loop(181) = {153, 150, 151, 152};
Line Loop(182) = {149, 146, 147, 148};
Line Loop(183) = {141, 135, 139, 140};
Line Loop(184) = {85, 82, 83, 84};
Line Loop(185) = {94, 98, 96, 97};
Line Loop(186) = {93, 90, 91, 92};
Line Loop(187) = {89, 86, 87, 88};
Line Loop(188) = {95, 99, 100, 101};
Line Loop(189) = {104, 105, 102, 103};
Line Loop(190) = {113, 110, 111, 112};
Line Loop(191) = {109, 106, 107, 108};
Plane Surface(192) = {156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191};
