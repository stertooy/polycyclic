gap> START_TEST( "IsMatrixRepresentation.tst" );

# The default images follow GeneratorsOfGroup, including redundant generators.
gap> A := AbelianPcpGroup( 1 );;
gap> M := [ [ 1, 1 ], [ 0, 1 ] ];;
gap> I := M^0;;
gap> H := Subgroup( A, [ A.1^2, A.1^4 ] );;
gap> IsMatrixRepresentation( H, [ M, M^2 ] );
true
gap> IsMatrixRepresentation( H, [ M, M ] );
false
gap> H := Subgroup( A, [ A.1, One( A ) ] );;
gap> IsMatrixRepresentation( H, [ M, I ] );
true
gap> IsMatrixRepresentation( H, [ M, M ] );
false

# Explicit generators must generate the specified source and belong to it.
gap> IsMatrixRepresentation( A, [ A.1 ], [ M ] );
true
gap> IsMatrixRepresentation( A, [ A.1, A.1^2 ], [ M, M^2 ] );
true
gap> IsMatrixRepresentation( A, [ A.1, A.1^2 ], [ M, M ] );
false
gap> IsMatrixRepresentation( A, [ A.1^2 ], [ M ] );
false
gap> H := Subgroup( A, [ A.1^2 ] );;
gap> IsMatrixRepresentation( H, [ A.1 ], [ M ] );
false
gap> IsMatrixRepresentation( A, [], [] );
false

# Empty lists and identity generators describe the trivial representation.
gap> T := AbelianPcpGroup( [] );;
gap> IsMatrixRepresentation( T, [] );
true
gap> IsMatrixRepresentation( T, [], [] );
true
gap> IsMatrixRepresentation( T, [ One( T ) ], [ I ] );
true
gap> IsMatrixRepresentation( T, [ One( T ) ], [ M ] );
false
gap> H := Subgroup( T, [ One( T ) ] );;
gap> IsMatrixRepresentation( H, [ I ] );
true
gap> IsMatrixRepresentation( H, [ M ] );
false

# Both forms reject unequal lengths, including extra matrices.
gap> IsMatrixRepresentation( A, [] );
Error, <gens> and <matrices> must be lists of same length
gap> IsMatrixRepresentation( A, [ M, I ] );
Error, <gens> and <matrices> must be lists of same length
gap> IsMatrixRepresentation( A, [ A.1, A.1^2 ], [ M ] );
Error, <gens> and <matrices> must be lists of same length
gap> IsMatrixRepresentation( T, [], [ I ] );
Error, <gens> and <matrices> must be lists of same length

# Explicit lists may reorder or omit dependent polycyclic generators.
gap> G := UnitriangularPcpGroup( 3, 0 );;
gap> x := [ [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 0, 1 ] ];;
gap> y := [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ] ];;
gap> z := Comm( x, y );;
gap> IsMatrixRepresentation( G, [ G.3, G.1, G.2 ], [ z, x, y ] );
true
gap> IsMatrixRepresentation( G, [ G.3, G.1, G.2 ], [ x, y, z ] );
false
gap> IsMatrixRepresentation( G, [ G.1, G.2 ], [ x, y ] );
true
gap> IsMatrixRepresentation( G, GeneratorsOfPcp( Pcp( G ) ), [ x, y, z ] );
true
gap> H := Subgroup( G, [ G.2, G.3 ] );;
gap> IsMatrixRepresentation( H, [ y, z ] );
true
gap> IsMatrixRepresentation( H, [ y, x ] );
false
gap> IsMatrixRepresentation( H, GeneratorsOfPcp( Pcp( H ) ), [ y, z ] );
true

# Quotients can retain an identity among their original group generators.
gap> N := SubgroupByIgs( G, [ G.3 ] );;
gap> F := ImagesSource( NaturalHomomorphismByNormalSubgroupNC( G, N ) );;
gap> GeneratorsOfGroup( F ) =
>        Concatenation( GeneratorsOfPcp( Pcp( F ) ), [ One( F ) ] );
true
gap> IsMatrixRepresentation( F, [ M, I, I ] );
true
gap> IsMatrixRepresentation( F, [ M, I, M ] );
false
gap> IsMatrixRepresentation( F, GeneratorsOfPcp( Pcp( F ) ), [ M, I ] );
true
gap> IsMatrixRepresentation( F, GeneratorsOfPcp( Pcp( F ) ),
>        [ M, TransposedMat( M ) ] );
false
gap> IsMatrixRepresentation( F, [ M, I ] );
Error, <gens> and <matrices> must be lists of same length

gap> STOP_TEST( "IsMatrixRepresentation.tst" );
