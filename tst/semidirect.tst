gap> START_TEST("Test of semidirect products");

#
gap> N := DihedralPcpGroup( 0 );;
gap> A := Group([
>    InnerAutomorphism( N, N.1 ),
>    InnerAutomorphism( N, N.2 ),
>    GroupHomomorphismByImagesNC( N, N, [ N.1, N.2 ], [ N.1*N.2, N.2^-1 ] )
> ]);;
gap> G := AbelianPcpGroup( [ 0 ] );;
gap> alpha := GroupHomomorphismByImagesNC( G, A, [ G.1 ], [ A.1*A.3 ] );;
gap> S := SemidirectProduct( G, alpha, N );
Pcp-group with orders [ 0, 2, 0 ]
gap> e1 := Embedding( S, 1 );;
gap> e2 := Embedding( S, 2 );;
gap> p := Projection( S );;
gap> e1 * p = IdentityMapping( G );
true
gap> ClosureGroup( Image( e1, G ), Image( e2, N ) ) = S;
true

# A finite permutation group acting nontrivially on an infinite pcp group.
# An action of order three distinguishes the action from its inverse.
gap> G := Group( (1,2,3) );;
gap> N := AbelianPcpGroup( [ 0, 0 ] );;
gap> a := GroupHomomorphismByImagesNC( N, N, [ N.1, N.2 ],
>                                    [ N.2, N.1^-1*N.2^-1 ] );;
gap> SetIsBijective( a, true );
gap> alpha := GroupHomomorphismByImagesNC( G, Group( [ a ] ), [ G.1 ], [ a ] );;
gap> S := SemidirectProduct( G, alpha, N );;
gap> IsPcpGroup( S );
true
gap> e1 := Embedding( S, 1 );;
gap> e2 := Embedding( S, 2 );;
gap> e1 * Projection( S ) = IdentityMapping( G );
true
gap> Image( e2, N.1 ) ^ Image( e1, G.1 ) = Image( e2, N.2 );
true

# The same action in the two-argument form, with N infinite.
gap> Aut := Group( [ a ] );;
gap> S := SemidirectProduct( Aut, N );;
gap> IsPcpGroup( S );
true
gap> Size( Aut );
3
gap> Embedding( S, 1 ) * Projection( S ) = IdentityMapping( Aut );
true
gap> Image( Embedding( S, 2 ), N.1 ) ^ Image( Embedding( S, 1 ), a )
>    = Image( Embedding( S, 2 ), N.2 );
true
gap> Kernel( Projection( S ) ) = Image( Embedding( S, 2 ) );
true

# A pcp group acting nontrivially on a finite permutation group.
gap> G := AbelianPcpGroup( [ 0 ] );;
gap> N := Group( (1,2,3,4,5) );;
gap> a := GroupHomomorphismByImages( N, N, [ N.1 ], [ N.1^2 ] );;
gap> alpha := GroupHomomorphismByImagesNC( G, Group( [ a ] ), [ G.1 ], [ a ] );;
gap> S := SemidirectProduct( G, alpha, N );;
gap> IsPcpGroup( S );
true
gap> e1 := Embedding( S, 1 );;
gap> e2 := Embedding( S, 2 );;
gap> e1 * Projection( S ) = IdentityMapping( G );
true
gap> Image( e2, N.1 ) ^ Image( e1, G.1 ) = Image( e2, N.1^2 );
true

# The two-argument form with N finite.
gap> N := AbelianPcpGroup( [ 5 ] );;
gap> a := GroupHomomorphismByImages( N, N, [ N.1 ], [ N.1^2 ] );;
gap> Aut := Group( [ a ] );;
gap> S := SemidirectProduct( Aut, N );;
gap> IsPcpGroup( S );
true
gap> Size( S );
20

#
gap> STOP_TEST( "semidirect.tst", 10000000);
