// Inputs
gridsize=0.04;
radius=0.1;
radius_int = 0.05;
a = 0.3;
phi = Asin(0.5/2);
alpha = Asin((Sin(phi)*radius)/radius_int);
len = 1;

xc = 0;
yc = 0;
zc = 0;

// Inner circle
Point(99995)={xc,yc,zc,gridsize};
Point(99996)={radius+xc,yc,zc,gridsize};
Point(99997)={-radius+xc,yc,zc,gridsize};
Point(99998)={xc,-radius+yc,zc,gridsize};
Point(99999)={xc,radius+yc,zc,gridsize};


// inner circle resonator
Point(999910)={radius_int+xc,yc,zc,gridsize};
Point(999911)={-radius_int+xc,yc,zc,gridsize};
Point(999912)={xc,-radius_int+yc,zc,gridsize};
Point(999913)={xc,radius_int+yc,zc,gridsize};

Rotate {{0, 0, 1}, {0, 0, 0}, phi} {
  Duplicata { Point{999911}; }
}
Rotate {{0, 0, 1}, {0, 0, 0}, -phi} {
  Duplicata { Point{999911}; }
}



Rotate {{0, 0, 1}, {xc, yc, zc}, phi} {
  Duplicata { Point{99997}; }
}
Rotate {{0, 0, 1}, {xc, yc, zc}, -phi} {
  Duplicata { Point{99997}; }
}


Point(999999)={-radius_int*Cos(alpha)+xc,radius_int*Sin(alpha)+yc,zc,gridsize};
Point(9999100)={-radius_int*Cos(alpha)+xc,-radius_int*Sin(alpha)+yc,zc,gridsize};

Circle(99995) = {99996, 99995, 999917};
Circle(99996) = {99996, 99995, 999916};
Circle(99997) = {999999, 99995, 999910};
Circle(99998) = {999910, 99995, 9999100};
Line(99999) = {999917, 999999};
Line(999910) = {999916, 9999100};


Line Loop(999912) = {99995, 99999, 99997, 99998, -999910, -99996};
Plane Surface(999913) = {-999912};
Delete {
  Point{99997, 999911, 999915, 999914, 99995};
}

Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{999913};
}
Extrude {0, 0, len} {
 Surface{999913};
}


Translate {a, 0, 0} {
  Duplicata { Surface{999924, 999932, 999928, 999945, 999940, 999936, 999944, 999913}; }
}
Translate {2*a, 0, 0} {
  Duplicata { Surface{999924, 999932, 999928, 999945, 999940, 999936, 999944, 999913}; }
}
Translate {3*a, 0, 0} {
  Duplicata { Surface{999924, 999932, 999928, 999945, 999940, 999936, 999944, 999913}; }
}
Translate {4*a, 0, 0} {
  Duplicata { Surface{999924, 999932, 999928, 999945, 999940, 999936, 999944, 999913}; }
}
Translate {0, a, 0} {
  Duplicata { Surface{999944, 999994, 999913, 999936, 999924, 999945, 999932, 999940, 999928, 999946, 999973, 999951, 999978, 999983, 999968, 999956, 999961, 999989, 1000011, 999984, 999999, 1000037, 1000075, 1000060, 1000087, 1000082, 1000049, 1000022, 1000027, 1000044, 1000016, 1000006, 1000021, 1000059, 1000054, 1000032, 1000070, 1000097, 1000092, 1000065}; }
}
Translate {0, 2*a, 0} {
  Duplicata { Surface{999944, 999994, 999913, 999936, 999924, 999945, 999932, 999940, 999928, 999946, 999973, 999951, 999978, 999983, 999968, 999956, 999961, 999989, 1000011, 999984, 999999, 1000037, 1000075, 1000060, 1000087, 1000082, 1000049, 1000022, 1000027, 1000044, 1000016, 1000006, 1000021, 1000059, 1000054, 1000032, 1000070, 1000097, 1000092, 1000065}; }
}
Translate {0, 3*a, 0} {
  Duplicata { Surface{999944, 999994, 999913, 999936, 999924, 999945, 999932, 999940, 999928, 999946, 999973, 999951, 999978, 999983, 999968, 999956, 999961, 999989, 1000011, 999984, 999999, 1000037, 1000075, 1000060, 1000087, 1000082, 1000049, 1000022, 1000027, 1000044, 1000016, 1000006, 1000021, 1000059, 1000054, 1000032, 1000070, 1000097, 1000092, 1000065}; }
}

Physical Surface(1000746) = {1000557, 1000341, 1000125, 999945, 1000616, 1000638, 1000645, 1000652, 999961, 999999, 1000037, 1000075, 1000220, 1000436, 1000213, 1000206, 1000184, 1000422, 1000400, 1000429};
Physical Surface(1000747) = {-1000540, -1000599, -1000324, -1000108, -1000167, -1000383, -1000488, -1000733, -1000517, -1000711, -1000704, -1000495, -999913, -999983, -1000021, -1000059, -1000097, -1000301, -1000279, -1000272};
Physical Surface(1000748) = {1000552, 1000564, 1000574, 1000547, 1000569, 1000530, 1000336, 1000348, 1000579, 1000589, 1000358, 1000331, 1000353, 1000611, 1000314, 1000584, 1000606, 1000120, 1000594, 1000132, 1000363, 1000373, 1000633, 1000142, 1000115, 1000623, 1000137, 1000395, 1000098, 1000368, 1000390, 999924, 1000535, 1000378, 1000628, 999932, 1000699, 1000147, 1000694, 1000157, 1000417, 999928, 999936, 1000407, 1000679, 999940, 1000179, 999944, 1000152, 1000684, 1000174, 1000319, 1000162, 1000412, 1000483, 999946, 1000723, 1000478, 1000674, 999951, 1000689, 1000201, 1000718, 1000191, 1000463, 999956, 999973, 1000468, 1000659, 999968, 1000103, 999978, 1000196, 1000745, 1000267, 1000507, 1000262, 1000458, 1000473, 999984, 1000728, 1000502, 1000664, 999989, 1000669, 1000247, 1000740, 1000252, 1000443, 999994, 1000011, 1000529, 1000006, 1000291, 1000016, 1000242, 1000257, 1000512, 1000286, 1000448, 1000453, 1000022, 1000524, 1000027, 1000227, 1000313, 1000032, 1000049, 1000044, 1000296, 1000054, 1000232, 1000237, 1000308, 1000060, 1000065, 1000070, 1000087, 1000082, 1000092};
