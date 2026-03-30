gap> START_TEST("Test for the centralizer algorithm");

#
gap> G := ExamplesOfSomePcpGroups( 5 );;
gap> H := Subgroup( G, [ G.1 ] );;
gap> C := Centralizer( H, G.2 );
Pcp-group with orders [ 0 ]
gap> C = Subgroup( G, [ G.4 ] );
true

#
gap> G := PcGroupToPcpGroup( PcGroupCode( 520, 16 ) );;
gap> g := G.1*G.3*G.4;;
gap> H := Subgroup( G,[ G.2, G.3, G.4 ] );;
gap> Centralizer( H, g );
Pcp-group with orders [ 2, 2 ]

#
gap> STOP_TEST( "cent.tst", 10000000);
