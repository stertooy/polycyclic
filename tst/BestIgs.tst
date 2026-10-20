gap> START_TEST("BestIgs.tst");

# BestIgs computes only an igs when no sequence has been stored.
gap> G := HeisenbergPcpGroup(1);;
gap> gens := [G.1^-2*G.2^8*G.3^20, G.2^-3*G.3^13, G.3^6];;
gap> U := Subgroup(G, gens);;
gap> [HasIgs(U), HasNgs(U), HasCgs(U)];
[ false, false, false ]
gap> igs := BestIgs(U);;
gap> [HasIgs(U), HasNgs(U), HasCgs(U), IsIdenticalObj(igs, Igs(U))];
[ true, false, false, true ]

# A later ngs or cgs takes precedence; the original Igs attribute is stable.
gap> ngs := Ngs(U);;
gap> [IsIdenticalObj(BestIgs(U), ngs), HasCgs(U)];
[ true, false ]
gap> cgs := Cgs(U);;
gap> [IsIdenticalObj(BestIgs(U), cgs), IsIdenticalObj(Igs(U), igs), igs = cgs];
[ true, true, false ]

# A stored cgs can also be used without computing an igs or ngs.
gap> U := Subgroup(G, gens);;
gap> cgs := Cgs(U);;
gap> [IsIdenticalObj(BestIgs(U), cgs), HasIgs(U), HasNgs(U)];
[ true, false, false ]

# The empty sequence is a valid result too.
gap> U := Subgroup(G, []);;
gap> BestIgs(U);
[  ]
gap> [IsIdenticalObj(BestIgs(U), Igs(U)), HasNgs(U), HasCgs(U)];
[ true, false, false ]

# Membership, closure, and both commutator branches agree for all choices.
# SubgroupByIgs preserves the nonnormed input sequence in these fixtures.
gap> check := function(G, gens)
> local U, V, samples, results, sequence;
> U := SubgroupByIgs(G, gens);
> V := Subgroup(G, [G.1]);
> samples := Concatenation([One(G), G.1, G.2], gens);
> results := [];
> for sequence in [Igs, Ngs, Cgs] do
>     sequence(U); sequence(V);
>     Add(results, [List(samples, x -> x in U),
>         Cgs(ClosureGroup(U, V)), Cgs(CommutatorSubgroup(U, V)),
>         Cgs(CommutatorSubgroup(U, U))]);
> od;
> return results[1] = results[2] and results[2] = results[3];
> end;;
gap> check(G, gens);
true

# Equal groups with identical stored igs can have different improved caches.
gap> U := SubgroupByIgs(G, gens);; V := SubgroupByIgs(G, gens);;
gap> Cgs(U);;
gap> D := CommutatorSubgroup(U, U);;
gap> [CommutatorSubgroup(U, V) = D, CommutatorSubgroup(V, U) = D];
[ true, true ]

gap> G := AbelianPcpGroup([12, 12]);;
gap> check(G, [G.1^5*G.2^2, G.2^6]);
true

# Centralizers, normalizers, and subgroup series also accept a cached cgs
# alongside a different stored igs. Compare independently constructed inputs.
gap> check := function(G, gens)
> local U, V, W, X, igs, before, series, fun, ser;
> U := SubgroupByIgs(G, gens);
> V := SubgroupByIgs(G, [gens[1]]);
> series := PcpSeries(U);
> before := [Cgs(Centralizer(U, V)), Cgs(Normalizer(U, V))];
> W := SubgroupByIgs(G, gens);
> X := SubgroupByIgs(G, [gens[1]]);
> igs := Igs(W);
> Cgs(W); Cgs(X);
> if igs = Cgs(W) or PcpSeries(W) <> series or
>    before <> [Cgs(Centralizer(W, X)), Cgs(Normalizer(W, X))] then
>     return false;
> fi;
> for fun in [EfaSeries, RefinedDerivedSeries, RefinedDerivedSeriesDown,
>             TorsionByPolyEFSeries] do
>     ser := fun(W);
>     if ser <> fail and not (ser[1] = W and IsTrivial(Last(ser)) and
>        ForAll([1..Length(ser)-1], i -> IsNormal(ser[i], ser[i+1]) and
>            IsAbelian(ser[i]/ser[i+1]))) then
>         return false;
>     fi;
> od;
> return IsIdenticalObj(Igs(W), igs);
> end;;
gap> G := DihedralPcpGroup(60);;
gap> check(G, [G.1*G.2^8, G.2^3]);
true
gap> G := DihedralPcpGroup(0);;
gap> check(G, [G.1*G.2^8, G.2^3]);
true
gap> G := HeisenbergPcpGroup(1);;
gap> check(G, [G.1^-2*G.2^8*G.3^20, G.2^-3*G.3^13, G.3^6]);
true

# Power series and torsion subgroups use the cached cgs as well.
gap> G := UnitriangularPcpGroup(3, 3);;
gap> U := SubgroupByIgs(G, [G.1*G.2^2, G.2, G.3]);;
gap> Cgs(U);; ser := PCentralSeries(U, 3);;
gap> List(ser, Size);
[ 27, 3, 1 ]
gap> PRump(U, 3) = ser[2];
true
gap> G := DirectProduct(HeisenbergPcpGroup(1), AbelianPcpGroup([4]));;
gap> U := SubgroupByIgs(G, [G.1*G.2^8, G.2, G.3, G.4]);;
gap> Cgs(U);; IsNilpotentGroup(U);;
gap> TorsionSubgroup(U) = Subgroup(G, [G.4]);
true

gap> STOP_TEST("BestIgs.tst", 1);
