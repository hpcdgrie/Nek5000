r=0.4;
Py=1.0;
Px=Py*3^0.5;
Point(1) = {Px/2.0,0,0};
Point(2) = {-Px/2.0,0,0};
Point(3) = {0,Py/2.0,0};
Point(4) = {0,-Py/2.0,0};
Point(5) = {Px/2.0-r,0,0};
Point(6) = {-Px/2.0+r,0,0};
Point(7) = {0,Py/2.0-r,0};
Point(8) = {0,-Py/2.0+r,0};
Point(9) = {(1.0/3.0)*Px/2.0,0,0};
Point(10) = {-(1.0/3.0)*Px/2.0,0,0};
Point(11) = {Px*0.5*(1.0-(r/Py)),Py*0.5*(r/Py),0};
Point(12) = {Px*0.5*(1.0-((Py-r)/Py)),Py*0.5*((Py-r)/Py),0};
Point(13) = {0.5*(Px*0.5*(1.0-(r/Py))+Px*0.5*(1.0-((Py-r)/Py))),0.5*(Py*0.5*(r/Py)+Py*0.5*((Py-r)/Py)),0};
Point(14) = {Px*0.5*(1.0-(r/Py)),-Py*0.5*(r/Py),0};
Point(15) = {Px*0.5*(1.0-((Py-r)/Py)),-Py*0.5*((Py-r)/Py),0};
Point(16) = {0.5*(Px*0.5*(1.0-(r/Py))+Px*0.5*(1.0-((Py-r)/Py))),-0.5*(Py*0.5*(r/Py)+Py*0.5*((Py-r)/Py)),0};
Point(17) = {-Px*0.5*(1.0-(r/Py)),Py*0.5*(r/Py),0};
Point(18) = {-Px*0.5*(1.0-((Py-r)/Py)),Py*0.5*((Py-r)/Py),0};
Point(19) = {-0.5*(Px*0.5*(1.0-(r/Py))+Px*0.5*(1.0-((Py-r)/Py))),0.5*(Py*0.5*(r/Py)+Py*0.5*((Py-r)/Py)),0};
Point(20) = {-Px*0.5*(1.0-(r/Py)),-Py*0.5*(r/Py),0};
Point(21) = {-Px*0.5*(1.0-((Py-r)/Py)),-Py*0.5*((Py-r)/Py),0};
Point(22) = {-0.5*(Px*0.5*(1.0-(r/Py))+Px*0.5*(1.0-((Py-r)/Py))),-0.5*(Py*0.5*(r/Py)+Py*0.5*((Py-r)/Py)),0};
Point(23) = {0,0,0};
Point(24) = {0.5*r,(1-0.5*r/((1.0/3.0)*Px/2.0))*Py*0.5,0};
Point(25) = {-0.5*r,(1-0.5*r/((1.0/3.0)*Px/2.0))*Py*0.5,0};
Point(26) = {0.5*r,-(1-0.5*r/((1.0/3.0)*Px/2.0))*Py*0.5,0};
Point(27) = {-0.5*r,-(1-0.5*r/((1.0/3.0)*Px/2.0))*Py*0.5,0};
//+
Circle(1) = {20, 2, 6};
//+
Circle(2) = {6, 2, 17};
//+
Circle(3) = {18, 3, 25};
//+
Circle(4) = {25, 3, 7};
//+
Circle(5) = {7, 3, 24};
//+
Circle(6) = {24, 3, 12};
//+
Circle(7) = {11, 1, 5};
//+
Circle(8) = {5, 1, 14};
//+
Circle(9) = {15, 4, 26};
//+
Circle(10) = {26, 4, 8};
//+
Circle(11) = {8, 4, 27};
//+
Circle(12) = {27, 4, 21};
//+
Line(13) = {17, 19};
//+
Line(14) = {18,19};
//+
Line(15) = {19, 10};
//+
Line(16) = {6, 10};
//+
Line(17) = {25,10};
//+
Line(18) = {20, 22};
//+
Line(19) = {22, 10};
//+
Line(20) = {21,22};
//+
Line(21) = {27, 10};
//+
Line(22) = {10, 23};
//+
Line(23) = {7, 23};
//+
Line(24) = {8,23};
//+
Line(25) = {23, 9};
//+
Line(26) = {24,9};
//+
Line(27) = {26,9};
//+
Line(28) = {9, 16};
//+
Line(29) = {15,16};
//+
Line(30) = {14,16};
//+
Line(31) = {5, 9};
//+
Line(32) = {13, 9};
//+
Line(33) = {12, 13};
//+
Line(34) = {11,13};
//+
Curve Loop(1) = {3, 17, -15, -14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {15, -16, 2, 13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, -19, -18, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {21, -19, -20, -12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {22, -23, -4, 17};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {22, -24, 11, 21};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {25, -26, -5, 23};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {25, -27, 10, 24};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {28, -29, 9, 27};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {8, 30, -28, -31};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {7, 31, -32, -34};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {33, 32, -26, 6};
//+
Plane Surface(12) = {12};
//+
Transfinite Curve {3, 2, 1, 15, 19, 4, 22, 11, 12, 10, 5, 25, 6, 7, 8, 9, 32, 28} = 6 Using Progression 1;
//+
Transfinite Curve {14, 13, 16, 17, 21, 20, 18, 23, 24, 26, 27, 31, 30, 29, 34, 33} = 15 Using Progression 1.2;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Transfinite Surface {7};
//+
Transfinite Surface {8};
//+
Transfinite Surface {9};
//+
Transfinite Surface {10};
//+
Transfinite Surface {11};
//+
Transfinite Surface {12};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Recombine Surface {5};
//+
Recombine Surface {6};
//+
Recombine Surface {7};
//+
Recombine Surface {8};
//+
Recombine Surface {9};
//+
Recombine Surface {10};
//+
Recombine Surface {11};
//+
Recombine Surface {12};
//+
Extrude {0, 0, 10} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{12}; Surface{11}; Surface{10}; Layers{50}; Recombine;
}
//+
Physical Surface("inlet", 1) = {1, 2, 3, 4, 6, 5, 7, 8, 12, 11, 10, 9};
//+
Physical Surface("outlet", 2) = {254, 276, 298, 232, 188, 210, 144, 166, 122, 100, 78, 56};
//+
Physical Surface("side1", 3) = {77, 55};
//+
Physical Surface("side2", 4) = {289, 223};
//+
Physical Surface("side3", 5) = {241, 275};
//+
Physical Surface("side4", 6) = {117, 95};
//+
Physical Surface("pin", 7) = {183, 253, 139, 43, 73, 99, 285, 121, 161, 205, 227, 263};
//+
Physical Volume("fluid", 8) = {4, 1, 6, 3, 2, 5, 7, 8, 9, 12, 11, 10};
