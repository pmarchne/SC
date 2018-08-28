// Inputs
gridsize=1/40;

//radius = DefineNumber[0.4, Name "R", Min 0.05, Max 0.45, Step 0.01];

DefineConstant[ radius = {0.45, Name "R", Min 0.05, Max 0.48, Step 0.01} ];


//radius=0.4;

a = 1;
xc=2;
yc=0.5;


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
Point(9)={radius+xc,yc,0,gridsize};
Point(7)={xc-radius,yc,0,gridsize};
Point(8)={xc,-radius+yc,0,gridsize};
Point(6)={xc,radius+yc,0,gridsize};

Circle(12) = {6, 5, 7};
Circle(13) = {7, 5, 8};
Circle(14) = {8, 5, 9};
Circle(15) = {9, 5, 6};


//Line Loop(28) = {4, 1, 2, 3};
//Line Loop(29) = {-12, -13, -14, -15};
Periodic Line {1}={3};
Periodic Line {2}={4};
//Plane Surface(30) = {28, 29};

Line Loop(32) = {2, 3, 4, 1};
Line Loop(33) = {12, 13, 14, 15};

Plane Surface(38) = {32, 33};
