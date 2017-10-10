lcar1 = .25;

Point(1) = {2,2,5,lcar1}; Point(2) = {-2,2,5,lcar1};
Point(3) = {2,-2,5,lcar1}; Point(4) = {-2,-2,5,lcar1};
Point(5) = {2,2,-5,lcar1}; Point(6) = {-2,2,-5,lcar1};
Point(7) = {2,-2,-5,lcar1}; Point(8) = {-2,-2,-5,lcar1};

// Instead of using included files, we now use a user-defined macro in order
// to carve some holes in the cube:

Macro CheeseHole

  // In the following commands we use the reserved variable name `newp', which
  // automatically selects a new point number. This number is chosen as the
  // highest current point number, plus one. (Note that, analogously to `newp',
  // the variables `newl', `news', `newv' and `newreg' select the highest number
  // amongst currently defined curves, surfaces, volumes and `any entities other
  // than points', respectively.)

  p1 = newp; Point(p1) = {x,  y,  z,  lcar3} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar3} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar3} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar3} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar3} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar3} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar3} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7}; c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4}; c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3}; c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6}; c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3}; c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6}; c12 = newreg; Circle(c12) = {p6,p1,p7};

  // We need non-plane surfaces to define the spherical holes. Here we use ruled
  // surfaces, which can have 3 or 4 sides:

  l1 = newreg; Line Loop(l1) = {c5,c10,c4};    Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};    Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1};  Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};   Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};   Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7};  Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12}; Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};   Surface(newreg) = {l8};

  // We then store the surface loops identification numbers in a list for later
  // reference (we will need these to define the final volume):

  theloops[t] = newreg ;

  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  thehole = newreg ;
  Volume(thehole) = theloops[t] ;

Return

t = 1; x = 0; y = 0; z = 0; r = 1; lcar3 = lcar1;

Call CheeseHole;

//+
Line(13) = {8, 7};
//+
Line(14) = {7, 5};
//+
Line(15) = {5, 6};
//+
Line(16) = {6, 8};
//+
Line(17) = {8, 4};
//+
Line(18) = {4, 3};
//+
Line(19) = {3, 1};
//+
Line(20) = {1, 2};
//+
Line(21) = {2, 4};
//+
Line(22) = {6, 2};
//+
Line(23) = {5, 1};
//+
Line(24) = {3, 7};
//+
Line Loop(28) = {20, -22, -15, 23};
//+
Plane Surface(29) = {28};
//+
Line Loop(29) = {16, 13, 14, 15};
//+
Plane Surface(30) = {29};
//+
Line Loop(30) = {17, 18, 24, -13};
//+
Plane Surface(31) = {30};
//+
Line Loop(31) = {19, 20, 21, 18};
//+
Plane Surface(32) = {31};
//+
Line Loop(32) = {22, 21, -17, -16};
//+
Plane Surface(33) = {32};
//+
Line Loop(33) = {14, 23, -19, 24};
//+
Plane Surface(34) = {33};

theloops[0] = newreg;
//+
Surface Loop(theloops[0]) = {32, 34, 30, 33, 29, 31};

thevolume = newreg;
//+
Volume(thevolume) = {theloops[]};

Physical Volume("flow_domain") = thevolume;
//+
Physical Surface("inlet") = {30};
//+
Physical Surface("outlet") = {32};
//+
Physical Surface("no_slip") = {31, 29, 34, 33, 26, 14, 22, 24, 28, 18, 20, 16};
