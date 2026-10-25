gap> START_TEST( "gh-00249.tst" );

#
# Fix matrix representations for pcp subgroups and reordered generators
# <https://github.com/gap-packages/polycyclic/issues/249>
#
gap> G := UnitriangularPcpGroup( 3, 0 );;
gap> H := Subgroup( G, [ G.1^2, G.2, G.3 ] );;
gap> IsNilpotentGroup( H );;
gap> phi := UnitriangularMatrixRepresentation( H );;
gap> Comm( H.1, H.2 )^phi = Comm( H.1^phi, H.2^phi );
true

# Example 1: respect relations of the subgroup.
gap> x := [ [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 0, 1 ] ];;
gap> y := [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ] ];;
gap> z := Comm( x, y );;
gap> IsMatrixRepresentation( H, [ x, y, z ] );
false
gap> IsMatrixRepresentation( H, [ x^2, y, z ] );
true

#
gap> G := AbelianPcpGroup( [ 2, 0 ] );;
gap> K := Subgroup( G, [ G.2, G.2^2 ] );;
gap> IsNilpotentGroup( K );;
gap> rep := UnitriangularMatrixRepresentation( K );;
gap> gensimgs := MappingGeneratorsImages( rep );;
gap> GroupHomomorphismByImages( K, Image(rep),
>        gensimgs[1], gensimgs[2] ) <> fail;
true

#
gap> G := AbelianPcpGroup( [ 2, 0 ] );;
gap> H := Subgroup( G, [ G.1 * G.2 ] );;
gap> IsNilpotentGroup( H );;
gap> RelativeOrdersOfPcp( Pcp( H ) );
[ 2, 0 ]
gap> phi := UnitriangularMatrixRepresentation( H );;
gap> Image( phi, H.1 ) ^ 2 = Image( phi, H.1 ^ 2 );
true

# Example 3: follow the input generator order.
gap> A := UnitriangularPcpGroup( 3, 0 );;
gap> G := Group( [ A.3, A.1, A.2 ] );;
gap> IsMatrixRepresentation( G, [ z, x, y ] );
true
gap> IsMatrixRepresentation( G, [ x, y, z ] );
false

#
gap> STOP_TEST( "gh-00249.tst" );
