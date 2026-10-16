gap> START_TEST("pcpelms.tst");

gap> G := AbelianPcpGroup([0,12,0,5,0]);;
gap> List(Concatenation([One(G)], GeneratorsOfGroup(G)), FactorOrder);
[ fail, 0, 12, 0, 5, 0 ]
gap> FactorOrder(G.2^3 * G.4^2);
12

gap> STOP_TEST("pcpelms.tst");
