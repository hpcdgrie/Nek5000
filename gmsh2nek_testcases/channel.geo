// Gmsh project created on Fri Sep 14 14:55:38 2018
//+
Point(1) = {-0.5, 0.5, 0, 1.0};
//+
Point(2) = {0.5, 0.5, 0, 1.0};
//+
Point(3) = {0.5, -0.5, 0, 1.0};
//+
Point(4) = {-0.5, -0.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {1, 4};
//+
Curve Loop(1) = {1, 2, 3, -4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 4, 3, 2} = 6 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

//+
Extrude {0, 0, 10} {
  Surface{1}; Layers{20}; Recombine;
}
//+
Transfinite Curve {12, 11, 20, 16} = 21 Using Progression 1;
//+
Transfinite Curve {6, 7, 9, 8} = 6 Using Progression 1;
//+
Transfinite Surface {13};
//+
Transfinite Surface {17};
//+
Transfinite Surface {21};
//+
Transfinite Surface {25};
//+
Transfinite Surface {26};
//+
Recombine Surface {13};
//+
Recombine Surface {25};
//+
Recombine Surface {13};
//+
Recombine Surface {26};
//+
Recombine Surface {25};
//+
Recombine Surface {25};
//+
Recombine Surface {26};
//+
Recombine Surface {21};
//+
Recombine Surface {17};
//+
Physical Surface("inlet", 1) = {1};
//+
Physical Surface("outlet", 2) = {26};
//+
Physical Surface("wall", 3) = {21, 25, 13, 17};
//+
Physical Volume("fluid", 4) = {1};
//+
Physical Surface("inlet", 1) += {1};
//+
Physical Surface("outlet", 2) += {26};
//+
Physical Surface("wall", 3) += {25, 21, 13, 17};
