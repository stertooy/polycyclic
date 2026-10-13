gap> START_TEST("Test for various former bugs");

#
gap> # The following used to trigger an error starting with:
gap> # "SolutionMat: matrix and vector incompatible called from"
gap> K:=AbelianPcpGroup([3,3,3]);;
gap> A:=Subgroup(K,[K.1]);;
gap> cr:=CRRecordBySubgroup(K,A);;
gap> ExtensionsCR(cr);;

#
# Comparing homomorphisms used to be broken
gap> K:=AbelianPcpGroup(1,[3]);;
gap> hom1:=GroupHomomorphismByImages(K,K,[K.1],[K.1]);;
gap> hom2:=GroupHomomorphismByImages(K,K,[K.1^2],[K.1^2]);;
gap> hom1=hom2;
true
gap> hom1=IdentityMapping(K);
true
gap> hom2=IdentityMapping(K);
true

#
gap> # The following incorrectly triggered an error at some point
gap> IsTorsionFree(ExamplesOfSomePcpGroups(5));
true

#
gap> # Verify IsGeneratorsOfMagmaWithInverses warnings are silenced
gap> IsGeneratorsOfMagmaWithInverses(GeneratorsOfGroup(ExamplesOfSomePcpGroups(5)));
true

#
gap> # Check for a bug reported 2012-01-19 by Robert Morse
gap> g := PcGroupToPcpGroup(PcGroupCode(61763227873,48));
Pcp-group with orders [ 2, 2, 2, 2, 3 ]
gap> # The next two commands used to trigger errors
gap> NonAbelianTensorSquare(Centre(g));
Pcp-group with orders [ 8 ]
gap> NonAbelianExteriorSquare(Centre(g));
Pcp-group with orders [  ]

#
gap> # Check for a bug reported 2012-01-19 by Robert Morse
gap> F := FreeGroup("x","y");
<free group on the generators [ x, y ]>
gap> x := F.1;; y := F.2;;
gap> G := F/[x^2/y^24, y^24, y^x/y^23];
<fp group on the generators [ x, y ]>
gap> iso := IsomorphismPcGroup(G);
[ x, y ] -> [ f1, f2*f5 ]
gap> iso1 := IsomorphismPcpGroup(Image(iso));
[ f1, f2, f3, f4, f5 ] -> [ g1, g2, g3, g4, g5 ]
gap> G := Image(iso*iso1);
Pcp-group with orders [ 2, 2, 2, 2, 3 ]
gap> # The next command used to trigger an error
gap> NonAbelianTensorSquare(Image(iso*iso1));
Pcp-group with orders [ 2, 2, 3, 2, 2, 2, 2 ]

#
gap> # The problem with the previous example is/was that Igs(G)
gap> # is set to a non-standard value. Experiment with that some more
gap> G := Parent(G);; # does nothing in GAP >= 4.13 -- possibly due to https://github.com/gap-system/gap/pull/5631
gap> igs := [ G.1, G.2*G.5, G.3*G.4*G.5^2, G.4*G.5, G.5 ];;
gap> G := Subgroup(G, igs);;
gap> SetIgs(G, igs);
gap> Igs(G);
[ g1, g2*g5, g3*g4*g5^2, g4*g5, g5 ]
gap> # Unfortunately, it seems that a lot of code that
gap> # really should be using Ngs or Cgs is using Igs incorrectly.
gap> # For example, direct products could return *invalid* embeddings:
gap> D := DirectProduct(G, G);
Pcp-group with orders [ 2, 2, 2, 2, 3, 2, 2, 2, 2, 3 ]
gap> hom:=Embedding(D,1);;
gap> mapi:=MappingGeneratorsImages(hom);;
gap> GroupHomomorphismByImages(Source(hom),Range(hom),mapi[1],mapi[2]) <> fail;
true
gap> hom:=Projection(D,1);;
gap> mapi:=MappingGeneratorsImages(hom);;
gap> GroupHomomorphismByImages(Source(hom),Range(hom),mapi[1],mapi[2]) <> fail;
true

#
gap> # Check for bug computing Schur extension of infinite cyclic groups,
gap> # found by Max Horn 2012-05-25
gap> G:=AbelianPcpGroup(1,[0]);
Pcp-group with orders [ 0 ]
gap> # The next command used to trigger an error
gap> SchurExtension(G);
Pcp-group with orders [ 0 ]

#
gap> # Check for bug computing Schur extensions of subgroups, found by MH 2012-05-25.
gap> G:=HeisenbergPcpGroup(2);
Pcp-group with orders [ 0, 0, 0, 0, 0 ]
gap> H:=Subgroup(G,[G.2^3*G.3^2, G.1^9]);
Pcp-group with orders [ 0, 0, 0 ]
gap> # The next command used to trigger an error
gap> SchurExtension(H);
Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ]

#
gap> # Check for bug computing Schur extensions of subgroups, found by MH 2012-05-25.
gap> G:=HeisenbergPcpGroup(2);
Pcp-group with orders [ 0, 0, 0, 0, 0 ]
gap> H:=Subgroup(G,[G.1, G.2]);
Pcp-group with orders [ 0, 0 ]
gap> # The next command used to trigger an error
gap> SchurExtension(H);
Pcp-group with orders [ 0, 0, 0 ]

#
gap> # Check for bug computing normalizer of two subgroups, found by MH 2012-05-30.
gap> # The problem was caused by incorrect resp. overly restrictive use of Parent().
gap> G:=HeisenbergPcpGroup(2);
Pcp-group with orders [ 0, 0, 0, 0, 0 ]
gap> A:=Subgroup(Subgroup(G,[G.2,G.3,G.4,G.5]), [G.3]);
Pcp-group with orders [ 0 ]
gap> B:=Subgroup(Subgroup(G,[G.1,G.4,G.5]), [G.4]);
Pcp-group with orders [ 0 ]
gap> Normalizer(A,B);
Pcp-group with orders [ 0 ]
gap> # The following used to trigger the error "arguments must have a common parent group"
gap> Normalizer(B,A);
Pcp-group with orders [ 0 ]

#
gap> # In polycyclic 2.9 and 2.10, the code for 2-cohomology computations was broken.
gap> G := UnitriangularPcpGroup(3,0);
Pcp-group with orders [ 0, 0, 0 ]
gap> mats := G!.mats;
[ [ [ 1, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 
  [ [ 1, 0, 0 ], [ 0, 1, 1 ], [ 0, 0, 1 ] ], 
  [ [ 1, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ] ]
gap> C := CRRecordByMats(G,mats);;
gap> cc := TwoCohomologyCR(C);;
gap> cc.factor.rels;
[ 2, 0, 0 ]
gap> c := cc.factor.prei[2];
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1 ]
gap> cc.gcb;
[ [ 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 ], 
  [ -1, 0, 1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1 ] ]
gap> cc.gcc;
[ [ 1, 0, 0, 0, 0, -2, -1, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 1, 0, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 1, 0, 0, -2, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 ] ]

#
gap> # LowerCentralSeriesOfGroup for non-nilpotent pcp-groups used to trigger
gap> # an infinite recursion
gap> G := PcGroupToPcpGroup(PcGroupCode(25,6));
Pcp-group with orders [ 2, 3 ]
gap> LowerCentralSeriesOfGroup(G);
[ Pcp-group with orders [ 2, 3 ], Pcp-group with orders [ 3 ] ]

#
# Test for a regression in TorsionSubgroup (reported by Sam Tertooy)
#
gap> TorsionSubgroup(AbelianPcpGroup([4,3]));
Pcp-group with orders [ 4, 3 ]

#
# PcpElementByExponents and PcpElementByGenExpList returned wrong elements
# for negative exponents on generators of finite relative order
#
gap> G := UnitriangularPcpGroup( 4, 2 );;
gap> c := Collector( G );;
gap> PcpElementByGenExpList( c, [ 5, 1, 1, 0 ] ) = G.5;
true
gap> PcpElementByGenExpList( c, [ 5, 1, 1, -1 ] ) = G.5 * G.1^-1;
true
gap> value := function( G, word )
>     local gens;
>     gens := GeneratorsOfGroup( G );
>     return Product( [ 1, 3 .. Length( word ) - 1 ],
>                     i -> gens[ word[i] ] ^ word[i+1] );
> end;;
gap> checkGenExp := function( G )
>     local c, n;
>     c := Collector( G );
>     n := NumberOfGenerators( c );
>     return ForAll( Tuples( [ 1 .. n ], 2 ), g ->
>         ForAll( Tuples( [ -3 .. 4 ], 2 ), e ->
>             PcpElementByGenExpList( c, [ g[1], e[1], g[2], e[2] ] )
>             = value( G, [ g[1], e[1], g[2], e[2] ] ) ) );
> end;;
gap> checkExp := function( G, range )
>     local c, n;
>     c := Collector( G );
>     n := NumberOfGenerators( c );
>     return ForAll( Tuples( range, n ), v -> PcpElementByExponents( c, v )
>         = value( G, Concatenation( List( [ 1 .. n ], i -> [ i, v[i] ] ) ) ) );
> end;;
gap> checkGenExp( G );  checkExp( G, [ -1 .. 2 ] );
true
true
gap> SetFilterObj( c, UseLibraryCollector );
gap> checkGenExp( G );  checkExp( G, [ -1 .. 2 ] );
true
true

# D16 = < a, b > with pcgs a, b, a^2, a^4: the power relation of the first
# generator skips the second one
gap> ftl := FromTheLeftCollector( 4 );;
gap> for i in [ 1 .. 4 ] do SetRelativeOrder( ftl, i, 2 ); od;
gap> SetPower( ftl, 1, [ 3, 1 ] );  SetPower( ftl, 3, [ 4, 1 ] );
gap> SetConjugate( ftl, 2, 1, [ 2, 1, 3, 1 ] );
gap> SetConjugate( ftl, 3, 2, [ 3, 1, 4, 1 ] );
gap> D := PcpGroupByCollector( ftl );;
gap> Exponents( PcpElementByExponents( ftl, [ 2, -1, 0, 0 ] ) );
[ 0, 1, 1, 1 ]
gap> checkGenExp( D );  checkExp( D, [ -5 .. 6 ] );
true
true
gap> SetFilterObj( ftl, UseLibraryCollector );
gap> checkGenExp( D );  checkExp( D, [ -5 .. 6 ] );
true
true

# generators of finite and of infinite relative order
gap> M := ExamplesOfSomePcpGroups( 5 );;
gap> checkGenExp( M );  checkExp( M, [ -2 .. 2 ] );
true
true
gap> SetFilterObj( Collector( M ), UseLibraryCollector );
gap> checkGenExp( M );  checkExp( M, [ -2 .. 2 ] );
true
true

#
# CollectWordOrFail requires non-zero exponents, positive on generators of
# finite relative order
#
gap> ev := [ 0, 0, 0, 0, 1, 0 ];;
gap> CollectWordOrFail( Collector( UnitriangularPcpGroup( 4, 2 ) ), ev, [ 1, -1 ] );
Error, Assertion failure
gap> CollectWordOrFail( Collector( UnitriangularPcpGroup( 4, 2 ) ), ev, [ 1, 0 ] );
Error, Assertion failure
gap> CollectWordOrFail( c, ev, [ 1, -1 ] );
Error, Assertion failure

#
# PcpGroupFpGroupPcPres rejected a consistent presentation with a redundant
# relator containing an inverse of a generator of finite relative order
#
gap> F := FreeGroup( 4 );;
gap> f := GeneratorsOfGroup( F );;
gap> rels := [ f[1]^2 / f[3], f[2]^2, f[3]^2 / f[4], f[4]^2,
>              f[2]^f[1] / ( f[2] * f[3] ), f[3]^f[1] / f[3], f[4]^f[1] / f[4],
>              f[3]^f[2] / ( f[3] * f[4] ), f[4]^f[2] / f[4], f[4]^f[3] / f[4],
>              f[2]^f[1] / f[2] / ( f[3] * f[4] ) ];;
gap> Q := PcpGroupFpGroupPcPres( F / rels );;
gap> Size( Q );
16

#
gap> STOP_TEST( "bugfix.tst" );
