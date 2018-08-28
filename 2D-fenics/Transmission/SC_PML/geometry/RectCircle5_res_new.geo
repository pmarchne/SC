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
Translate {0, Ty, 0} {
  Duplicata { Line{77, 202, 201, 200, 199, 198, 197, 196, 88, 203, 66, 55, 12, 11, 211, 219, 218, 217, 216, 215, 214, 213, 212, 210, 209, 208, 207, 206, 205, 204}; }
}
Translate {0, -Ty, 0} {
  Duplicata { Line{77, 202, 201, 200, 199, 198, 197, 196, 88, 203, 66, 55, 12, 11, 211, 219, 218, 217, 216, 215, 214, 213, 212, 210, 209, 208, 207, 206, 205, 204}; }
}
Translate {0, 2*Ty, 0} {
  Duplicata { Line{77, 202, 201, 200, 199, 198, 197, 196, 88, 203, 66, 55, 12, 11, 211, 219, 218, 217, 216, 215, 214, 213, 212, 210, 209, 208, 207, 206, 205, 204}; }
}
Translate {0, -2*Ty, 0} {
  Duplicata { Line{77, 202, 201, 200, 199, 198, 197, 196, 88, 203, 66, 55, 12, 11, 211, 219, 218, 217, 216, 215, 214, 213, 212, 210, 209, 208, 207, 206, 205, 204}; }
}Line Loop(340) = {1, 2, 3, 4};
Line Loop(341) = {226, 225, 227, 222, -224, -223};
Line Loop(342) = {236, 237, -235, -240, -238, -239};
Line Loop(343) = {234, -241, -245, -243, -244, 242};
Line Loop(344) = {233, 220, 228, -232, -230, 231};
Line Loop(345) = {229, 249, 221, 246, -248, -247};
Line Loop(346) = {219, -217, -218, 215, 216, 214};
Line Loop(347) = {257, 252, -254, -253, 256, 255};
Line Loop(348) = {251, 276, -278, -277, 259, 279};
Line Loop(349) = {250, 258, -262, -260, 261, 263};
Line Loop(350) = {88, -12, -66, 55, 11, 77};
Line Loop(351) = {196, 201, -199, -200, 197, 198};
Line Loop(352) = {206, 205, -207, -202, -204, -203};
Line Loop(353) = {213, -211, -212, 209, 210, 208};
Line Loop(354) = {309, 281, 306, -308, -307, 289};
Line Loop(355) = {319, 339, 311, 336, -338, -337};
Line Loop(356) = {318, -322, -320, 321, 323, 310};
Line Loop(357) = {317, 312, -314, -313, 316, 315};
Line Loop(358) = {305, 301, -294, -302, 304, 303};
Line Loop(359) = {300, 295, -297, -296, 299, 298};
Line Loop(360) = {330, 325, -327, -326, 329, 328};
Line Loop(361) = {335, 331, -324, -332, 334, 333};
Line Loop(362) = {268, 270, 265, -267, -266, 269};
Line Loop(363) = {275, 271, -264, -272, 274, 273};
Line Loop(364) = {293, 280, 288, -292, -290, 291};
Line Loop(365) = {287, 282, -284, -283, 286, 285};
Plane Surface(366) = {340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365};
