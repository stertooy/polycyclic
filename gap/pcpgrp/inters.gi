#############################################################################
##
#W  grpint.gi                  Polycyc                           Bettina Eick
##

#############################################################################
##
#F NormalIntersection( N, U ) . . . . . . . . . . . . . . . . . . .  U \cap N
##
## The whole idea here is that the intersection U \cap N equals the kernel of
## the natural homomorphism \phi : U \to UN/N ; see also section 8.8.1 of
## the "Handbook of computational group theory".
InstallMethod( NormalIntersection, "for pcp groups",
               IsIdenticalObj, [IsPcpGroup, IsPcpGroup],
function( N, U )
    local G, igs, igsN, igsU, n, s, I, UN, q, Q, gens, imgs, phi;

	# get common overgroup of N and U
	G := PcpGroupByCollector( Collector( N ) );

    igs  := Igs(G);
    igsN := Cgs( N );
    igsU := Cgs( U );
    n    := Length( igs );

    # if N or U is trivial
    if Length( igsN ) = 0 then
        return N;
    elif Length( igsU ) = 0 then
        return U;
    fi;

    # if N or U are equal to G
    if Length( igsN ) = n and ForAll(igsN, x -> LeadingExponent(x) = 1) then
        return U;
    elif Length(igsU) = n and ForAll(igsU, x -> LeadingExponent(x) = 1) then
        return N;
    fi;

    # if N is a tail, we can read off the result directly
    s := Depth( igsN[1] );
    if Length( igsN ) = n-s+1 and
       ForAll( igsN, x -> LeadingExponent(x) = 1 ) then
        I := Filtered( igsU, x -> Depth(x) >= s );
        return SubgroupByIgs( G, I );
    fi;

    # otherwise compute
    UN := ClosureGroup( U, N );
    q := NaturalHomomorphismByNormalSubgroupNC( UN, N );
    Q := Range( q );
    gens := GeneratorsOfGroup( U );
    imgs := List( gens, u -> u^q );
    phi := GroupHomomorphismByImagesNC( U, Q, gens, imgs );
    return KernelOfMultiplicativeGeneralMapping( phi );
end );

#############################################################################
##
#M Intersection( N, U )
##
InstallMethod( Intersection2, "for pcp groups",
               IsIdenticalObj, [IsPcpGroup, IsPcpGroup],
function( U, V )
    # check for trivial cases
    if IsInt(Size(U)) and IsInt(Size(V)) then
        if IsInt(Size(V)/Size(U)) and ForAll(Igs(U), x -> x in V ) then
            return U;
        elif Size(V)<Size(U) and IsInt(Size(U)/Size(V))
             and ForAll( Igs(V), x -> x in U ) then
            return V;
        fi;
    fi;

    # test if one the groups is known to be normal
    if IsNormal( V, U ) then
        return NormalIntersection( U, V );
    elif IsNormal( U, V ) then
        return NormalIntersection( V, U );
    fi;

    if IsFinite( U ) or IsFinite( V ) then
        TryNextMethod();
    fi;
    Error("sorry: intersection for non-normal groups not yet installed");
end );
