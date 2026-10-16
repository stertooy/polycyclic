gap> START_TEST("Test of conjugacy");

#
gap> G := ExamplesOfSomePcpGroups(5);;
gap> H := Subgroup( G, [ G.1 * G.2 ^ -3 * G.3 ^ 2 * G.4 ^ -1,
> G.1 * G.2 ^ -1 * G.3 ^ -1 * G.4 ] );;
gap> g := G.2 ^ 2 * G.3 ^ -1 * G.4 ^ -1;;
gap> h := G.2 ^ 2 * G.3 ^ -1 * G.4;;
gap> r := RepresentativeAction( H, g, h );
fail
gap> IsConjugate( H, h, g );
false
gap> r := RepresentativeAction( G, g, h );;
gap> g ^ r = h;
true
gap> IsConjugate( G, h, g );
true

# Subgroup conjugacy in a finite overgroup inside an infinite pcp group.
gap> P := DirectProduct(DihedralPcpGroup(6), AbelianPcpGroup(1));;
gap> IsFinite(P);
false
gap> H := Subgroup(P, [P.1, P.2]);;
gap> G := Subgroup(P, [P.2]);;
gap> U := Subgroup(P, [P.1]);;
gap> V := U^P.2;;
gap> r := RepresentativeAction(H, U, V);;
gap> r in H and U^r = V;
true
gap> RepresentativeAction(H, U, G);
fail
gap> r := RepresentativeAction(G, U, V);;
gap> r in G and U^r = V;
true
gap> IsConjugate(G, U, V);
true
gap> RepresentativeAction(U, U, V);
fail
gap> IsConjugate(U, U, V);
false

# Conjugating lines in an elementary abelian layer, with neither subgroup
# contained in the conjugating group.
gap> A := AlternatingGroup(4);;
gap> iso := IsomorphismPcpGroup(A);;
gap> P := Range(iso);;
gap> U := Subgroup(P, [Image(iso, (1,2)(3,4))]);;
gap> V := Subgroup(P, [Image(iso, (1,3)(2,4))]);;
gap> C := Subgroup(P, [Image(iso, (1,2,3))]);;
gap> r := RepresentativeAction(C, U, V);;
gap> r in C and U^r = V and IsConjugate(C, U, V);
true
gap> RepresentativeAction(U, U, V);
fail

# Subgroup conjugacy in an infinite overgroup
gap> P := DihedralPcpGroup(0);;
gap> U := Subgroup(P, [P.1]);;
gap> V := U^P.2;;
gap> Size(ClosureGroup(U, V));
infinity
gap> RepresentativeAction(U, U, V);
fail
gap> IsConjugate(P, U, V);
true
gap> RepresentativeAction(P, U, V);
g2
gap> RepresentativeAction(P, U, Subgroup(P, [P.1*P.2^3]));
fail
gap> RepresentativeAction(P, U, Subgroup(P, [P.2]));
fail

# Subgroup conjugacy through free abelian layers
gap> G := ExamplesOfSomePcpGroups(5);;
gap> U := Subgroup(G, [G.1*G.3^-1, G.2^2*G.4]);;
gap> g := G.1^2*G.2^-1*G.3*G.4^-2;;
gap> V := U^g;;
gap> r := RepresentativeAction(G, U, V);;
gap> U^r = V;
true
gap> RepresentativeAction(G, U, Subgroup(G, [G.1, G.2^2*G.4]));
fail
gap> H := Subgroup(G, [G.1*G.2, G.3^2, G.4^3]);;
gap> V := U^(G.1*G.2*G.3^2);;
gap> r := RepresentativeAction(H, U, V);;
gap> r in H and U^r = V;
true
gap> IsConjugate(H, U, V);
true

# Trivial group
gap> T := TrivialSubgroup(DihedralPcpGroup(0));;
gap> RepresentativeAction(T, T, T);
id

# Lattices in a free abelian layer swapped by the action.
gap> coll := FromTheLeftCollector(3);;
gap> SetConjugateNC(coll, 2, 1, [3, 1]);;
gap> SetConjugateNC(coll, 3, 1, [2, 1]);;
gap> UpdatePolycyclicCollector(coll);;
gap> G := PcpGroupByCollectorNC(coll);;
gap> U := Subgroup(G, [G.2^2, G.3]);;
gap> V := Subgroup(G, [G.2, G.3^2]);;
gap> RepresentativeAction(G, U, V);
g1
gap> RepresentativeAction(Subgroup(G, [G.1]), U, V);
g1
gap> RepresentativeAction(Subgroup(G, [G.1^2]), U, V);
fail
gap> RepresentativeAction(G, U, Subgroup(G, [G.2*G.3, G.3^2]));
fail
gap> RepresentativeAction(G, Subgroup(G, [G.2^6, G.3^3]), Subgroup(G, [G.2^3, G.3^6]));
g1

# The stabilizer becomes trivial within the module series of a layer
gap> P := BurdeGrunewaldPcpGroup(1,2);;
gap> U := Subgroup(P, [P.2^-1*P.4*P.7*P.10*P.11^-1,
> P.3^-1*P.4^-1*P.5^-1*P.6*P.9]);;
gap> g := P.1^-1*P.2*P.3*P.4^-1*P.5^-1*P.6*P.7^-1*P.8*P.10^-1;;
gap> r := RepresentativeAction(P, U, U^g);;
gap> U^r = U^g;
true

# A restricted conjugating group need not contain either finite subgroup.
gap> P := DihedralPcpGroup(0);;
gap> C := Subgroup(P, [P.2^2]);;
gap> U := Subgroup(P, [P.1]);;
gap> V := U^(P.2^2);;
gap> r := RepresentativeAction(C, U, V);;
gap> r in C and U^r = V and IsConjugate(C, U, V);
true
gap> RepresentativeAction(C, U, U^P.2);
fail
gap> IsConjugate(C, U, U^P.2);
false

# Integral quotient lifts must stay in the conjugating subgroup.
gap> P := HeisenbergPcpGroup(1);;
gap> C := Subgroup(P, [P.1^2*P.3^3]);;
gap> U := Subgroup(P, [P.2]);;
gap> V := U^(P.1^2);;
gap> P.1^2 in C;
false
gap> r := RepresentativeAction(C, U, V);;
gap> r in C and U^r = V and IsConjugate(C, U, V);
true
gap> RepresentativeAction(C, U, U^P.1);
fail
gap> IsConjugate(C, U, U^P.1);
false

# The free abelian layer modulo the intersection has order 12. Its
# successive prime factors distinguish the two complement parity classes.
gap> G := ExamplesOfSomePcpGroups(3);;
gap> U := Subgroup(G, [G.1, G.2^12]);;
gap> V := U^G.2;;
gap> W := Subgroup(G, [G.1*G.2, G.2^12]);;
gap> r := RepresentativeAction(G, U, V);;
gap> U^r = V and IsConjugate(G, U, V);
true
gap> RepresentativeAction(G, U, W);
fail
gap> C := Subgroup(G, [G.1, G.2^2]);;
gap> IsConjugate(C, U, V);
false
gap> IsConjugate(C, U, U^(G.2^2));
true

# Mixed first cohomology C3 x Z requires transporting the torsion class
# after transporting the free part, including negative acting exponents.
gap> coll := FromTheLeftCollector(6);;
gap> SetConjugateNC(coll, 3, 1, [3, 1, 4, -1, 6, 1]);;
gap> SetConjugateNC(coll, 3, 2, [3, 1, 6, 1]);;
gap> SetConjugateNC(coll, 4, 2, [4, -1]);;
gap> SetConjugateNC(coll, 5, 2, [5, -1]);;
gap> SetConjugateNC(coll, 4, 3, [5, 1]);;
gap> SetConjugateNC(coll, 5, 3, [4, -1, 5, -1]);;
gap> SetConjugateNC(coll, 6, 2, [4, -2, 6, 1]);;
gap> SetConjugateNC(coll, 6, 3, [4, -1, 5, 1, 6, 1]);;
gap> UpdatePolycyclicCollector(coll);;
gap> G := PcpGroupByCollectorNC(coll);;
gap> U := Subgroup(G, [G.3]);;
gap> ForAll([G.1, G.2, G.1^-1, G.2^-1, G.1*G.2^-1], function(g)
> local V, r;
> V := U^g;
> r := RepresentativeAction(G, U, V);
> return r <> fail and r in G and U^r = V and IsConjugate(G, U, V);
> end);
true

gap> W := Subgroup(G, [G.3*G.4^2]);;
gap> RepresentativeAction(G, U, W);
fail
gap> IsConjugate(G, U, W);
false

# Mixed torsion/free cocycle actions must handle inverse words both in
# the transporter and in generators stabilizing the free part.
gap> coll := FromTheLeftCollector(6);;
gap> SetConjugate(coll, 3, 1, [3, 1, 4, -1, 6, 1]);;
gap> SetConjugate(coll, 3, 2, [3, 1, 4, -1, 6, 1]);;
gap> SetConjugate(coll, 4, 1, [5, 1]);;
gap> SetConjugate(coll, 5, 1, [4, -1, 5, 5]);;
gap> SetConjugate(coll, 6, 1, [4, -1, 5, 1, 6, 1]);;
gap> SetConjugate(coll, 4, 2, [4, -1, 5, 5]);;
gap> SetConjugate(coll, 5, 2, [4, -5, 5, 24]);;
gap> SetConjugate(coll, 6, 2, [4, -2, 5, 5, 6, 1]);;
gap> SetConjugate(coll, 4, 3, [4, -115, 5, 551]);;
gap> SetConjugate(coll, 5, 3, [4, -551, 5, 2640]);;
gap> SetConjugate(coll, 6, 3, [4, -116, 5, 551, 6, 1]);;
gap> G := PcpGroupByCollector(coll);;
gap> H := Subgroup(G, [G.1, G.3, G.4, G.5, G.6]);;
gap> U := Subgroup(G, [G.3]);;
gap> V := U^(G.1^-1);;
gap> r := RepresentativeAction(H, U, V);;
gap> r <> fail and r in H and U^r = V and IsConjugate(H, U, V);
true
gap> U := Subgroup(G, [G.3*G.4]);;
gap> V := U^(G.1*G.2^-1);;
gap> r := RepresentativeAction(G, U, V);;
gap> r <> fail and r in G and U^r = V and IsConjugate(G, U, V);
true

#
gap> STOP_TEST( "conj.tst", 10000000);
