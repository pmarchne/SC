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
Translate {0, 4*Ty, 0} {
  Duplicata { Line{16, 21, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}; }
}
Translate {0, -4*Ty, 0} {
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
Line Loop(194) = {5, 6, 7, 8, 9, -3, -2, -1};
Plane Surface(195) = {194};
Line Loop(196) = {4, 1, 2, 3};
Line Loop(197) = {48, 49, 46, 47};
Line Loop(198) = {57, 54, 58, 56};
Line Loop(199) = {55, 59, 60, 61};
Line Loop(200) = {53, 50, 51, 52};
Line Loop(201) = {45, 42, 43, 44};
Line Loop(202) = {41, 35, 39, 40};
Line Loop(203) = {67, 68, 69, 66};
Line Loop(204) = {76, 77, 74, 78};
Line Loop(205) = {75, 79, 80, 81};
Line Loop(206) = {73, 70, 71, 72};
Line Loop(207) = {38, 36, 37, 34};
Line Loop(208) = {65, 62, 63, 64};
Line Loop(209) = {18, 19, 20, 21};
Line Loop(210) = {17, 14, 15, 16};
Line Loop(211) = {28, 29, 26, 27};
Line Loop(212) = {33, 30, 31, 32};
Line Loop(213) = {25, 22, 23, 24};
Line Loop(214) = {125, 122, 123, 124};
Line Loop(215) = {154, 158, 156, 157};
Line Loop(216) = {133, 130, 131, 132};
Line Loop(217) = {129, 126, 127, 128};
Line Loop(218) = {155, 159, 160, 161};
Line Loop(219) = {121, 115, 119, 120};
Line Loop(220) = {118, 116, 117, 114};
Line Loop(221) = {164, 165, 162, 163};
Line Loop(222) = {173, 170, 171, 172};
Line Loop(223) = {169, 166, 167, 168};
Line Loop(224) = {86, 87, 88, 89};
Line Loop(225) = {95, 99, 100, 101};
Line Loop(226) = {94, 98, 96, 97};
Line Loop(227) = {93, 90, 91, 92};
Line Loop(228) = {85, 82, 83, 84};
Line Loop(229) = {105, 102, 103, 104};
Line Loop(230) = {113, 110, 111, 112};
Line Loop(231) = {109, 106, 107, 108};
Line Loop(232) = {135, 139, 140, 141};
Line Loop(233) = {147, 148, 149, 146};
Line Loop(234) = {145, 142, 143, 144};
Line Loop(235) = {138, 136, 137, 134};
Line Loop(236) = {181, 175, 179, 180};
Line Loop(237) = {193, 190, 191, 192};
Line Loop(238) = {189, 186, 187, 188};
Line Loop(239) = {185, 182, 183, 184};
Line Loop(240) = {178, 176, 177, 174};
Line Loop(241) = {153, 150, 151, 152};
Plane Surface(242) = {196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241};
