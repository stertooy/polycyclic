gap> START_TEST( "fallback.tst" );
gap> useAlnuth := POLYCYCLIC_USE_ALNUTH;; POLYCYCLIC_USE_ALNUTH := false;;

# The first batch of primes finds no relations, so verification must retry.
gap> G := AbelianPcpGroup( [0,0] );;
gap> a := [[2,1],[1,1]];;
gap> StabilizerIntegralAction( G, [a^12,a^8], [1,0] ) =
> Subgroup( G, [G.1^2*G.2^-3] );
true

# Independent units in a cubic field: the bound must handle degree three
# and a primitive polynomial with constant coefficient -1.
gap> x := Indeterminate( Rationals );;
gap> a := CompanionMat( x^3-3*x+1 );;
gap> IsTrivial( StabilizerIntegralAction( G, [a,a-a^0], [1,0,0] ) );
true

# Neither generator is primitive in this biquadratic field. A primitive
# element must be constructed from the algebra basis and need not be a unit.
gap> a := [[2,1],[1,1]]^4;; b := [[3,2],[4,3]]^4;;
gap> mats := [KroneckerProduct(a,IdentityMat(2)),
> KroneckerProduct(IdentityMat(2),b)];;
gap> IsTrivial( StabilizerIntegralAction( G, mats, [1,0,0,0] ) );
true

gap> POLYCYCLIC_USE_ALNUTH := useAlnuth;;
gap> STOP_TEST( "fallback.tst" );
