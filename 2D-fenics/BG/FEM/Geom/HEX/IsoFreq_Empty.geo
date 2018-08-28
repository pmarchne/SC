gridsize=0.5;
xmax = 4/3*Pi;
ymax = 3.62759873;
xmid = 2.0943951; 

Point(1)={xmid,-ymax,0,gridsize};
Point(2)={xmax,0,0,gridsize};
Point(3)={xmid,ymax,0,gridsize};
Point(4)={-xmid,ymax,0,gridsize};
Point(5)={-xmax,0,0,gridsize};
Point(6)={-xmid,-ymax,0,gridsize};

Line(1) = {6, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};

Line Loop(7) = {3, 4, 5, 6, 1, 2};
Plane Surface(8) = {7};