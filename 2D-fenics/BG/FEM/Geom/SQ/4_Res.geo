h = 0.025;
Point(1) = {0.5, -0.5, 0, h};
Point(2) = {0.5, 0.5, 0, h};
Point(3) = {-0.5, 0.5, 0, h};
Point(4) = {-0.5, -0.5, 0, h};
Point(6) = {-0.2599999999999998, 0.4579898987322333, 0, h};
Point(7) = {-0.06201010126776652, 0.26, 0, h};
Point(8) = {-0.2599999999999998, 0.06201010126776668, 0, h};
Point(9) = {-0.4579898987322331, 0.26, 0, h};
Point(60) = {-0.2599999999999998, 0.4014213562373095, 0, h};
Point(70) = {-0.09029437251522854, 0.2317157287525381, 0, h};
Point(80) = {-0.2317157287525378, 0.09029437251522864, 0, h};
Point(90) = {-0.401421356237309, 0.26, 0, h};
Point(166) = {0.2599999999999998, 0.06201010126776668, 0, h};
Point(167) = {0.4579898987322331, 0.26, 0, h};
Point(170) = {0.401421356237309, 0.26, 0, h};
Point(171) = {0.2317157287525378, 0.09029437251522864, 0, h};
Point(174) = {0.2599999999999998, 0.4014213562373095, 0, h};
Point(179) = {0.2599999999999998, 0.4579898987322333, 0, h};
Point(182) = {0.09029437251522854, 0.2317157287525381, 0, h};
Point(187) = {0.06201010126776652, 0.26, 0, h};
Point(188) = {-0.2599999999999998, -0.06201010126776668, 0, h};
Point(189) = {-0.4579898987322331, -0.26, 0, h};
Point(192) = {-0.401421356237309, -0.26, 0, h};
Point(193) = {-0.2317157287525378, -0.09029437251522864, 0, h};
Point(196) = {-0.2599999999999998, -0.4014213562373095, 0, h};
Point(201) = {-0.2599999999999998, -0.4579898987322333, 0, h};
Point(204) = {-0.09029437251522854, -0.2317157287525381, 0, h};
Point(209) = {-0.06201010126776652, -0.26, 0, h};
Point(210) = {0.09029437251522854, -0.2317157287525381, 0, h};
Point(211) = {0.2599999999999998, -0.4014213562373095, 0, h};
Point(214) = {0.2599999999999998, -0.4579898987322333, 0, h};
Point(215) = {0.06201010126776652, -0.26, 0, h};
Point(222) = {0.2317157287525378, -0.09029437251522864, 0, h};
Point(223) = {0.2599999999999998, -0.06201010126776668, 0, h};
Point(226) = {0.401421356237309, -0.26, 0, h};
Point(231) = {0.4579898987322331, -0.26, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(10) = {6, 7};
Line(11) = {7, 70};
Line(12) = {70, 60};
Line(13) = {60, 90};
Line(14) = {90, 80};
Line(15) = {80, 8};
Line(16) = {8, 9};
Line(17) = {9, 6};
Line(40) = {166, 167};
Line(41) = {170, 171};
Line(42) = {174, 170};
Line(43) = {167, 179};
Line(44) = {182, 174};
Line(45) = {179, 187};
Line(46) = {187, 182};
Line(47) = {171, 166};
Line(48) = {188, 189};
Line(49) = {192, 193};
Line(50) = {196, 192};
Line(51) = {189, 201};
Line(52) = {204, 196};
Line(53) = {201, 209};
Line(54) = {209, 204};
Line(55) = {193, 188};
Line(56) = {210, 211};
Line(57) = {214, 215};
Line(58) = {215, 210};
Line(59) = {222, 223};
Line(60) = {226, 222};
Line(61) = {223, 231};
Line(62) = {211, 226};
Line(63) = {231, 214};
Line Loop(6) = {4,1,2,3,-16,-15,-14,-13,-12,-11,-10,-17,-40,-41,-42,-43,-44,-45,-46,-47,-48,-49,-50,-51,-52,-53,-54,-55,-56,-57,-58,-59,-60,-61,-62,-63};
Plane Surface(6) = {6};