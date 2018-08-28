// Inputs

gridsize=1/13;
a=10;
b=10;
r=0.0;
L_pml=8;



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


Line(1) = {7, 993};
Line(2) = {993, 992};
Line(3) = {992, 8};
Line(4) = {8, 7};
Line(5) = {7, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 4};
Line(9) = {4, 8};
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(11) = {10};
Line Loop(12) = {5, 6, 7, 8, 9, -3, -2, -1};
Plane Surface(13) = {12};
