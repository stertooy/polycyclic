#############################################################################
##
#W  pcppara.gi                   Polycyc                         Bettina Eick
#W                                                              Werner Nickel
##
## Parallel versions of the non-commuatative gauss algorithm.
##

#############################################################################
##
#F UpdateCounterPara( ind, c )  . . . . . . . . . . . . small help function
##
BindGlobal( "UpdateCounterPara", function( ind, c )
    local i;
    i := c - 1;
    while i > 0 and not IsBool(ind[i]) and LeadingExponent(ind[i]) = 1 do
        i := i - 1;
    od;
    return i + 1;
end );


#############################################################################
##
#F AddToIgsParallel( <pcs>, <gens>, <ppcs>, <pgens> )
##
## This function adds the elements in <gens> to the induced pcs <pcs>.
## It acts simultaneously on <pcs> and <ppcs> as well as <gens> and <pgens>.
##
InstallGlobalFunction( AddToIgsParallel,
function( pcs, gens, ppcs, pgens )
    local coll, rels, n, todo, tododo, ind, indd, g, gg, d, h, hh, k,
          eg, eh, e, changed, c, i, r, sub, val, j, f, a, b;

    if Length( gens ) = 0 then return [pcs, ppcs]; fi;

    # get information
    coll := Collector( gens[1] );
    rels := RelativeOrders( coll );
    n    := NumberOfGenerators( coll );
    c    := n+1;

    # set up
    ind  := ListWithIdenticalEntries(n, false);
    indd := ListWithIdenticalEntries(n, false);
    for i in [1..Length(pcs)] do
        d := Depth( pcs[i] );
        ind[d]  := pcs[i];
        indd[d] := ppcs[i];
    od;

    # do a reduction step
    c      := TailLimit(ind, c);
    sub    := Filtered([1..Length(gens)], i -> Depth(gens[i]) < c);
    sub    := Filtered(sub, i -> not gens[i] in gens{[1..i-1]}); # Essentially Set(todo)
    todo   := gens{sub};
    tododo := pgens{sub};
    val    := List(todo, x -> IGSValFun(x));

    # loop over to-do list until it is empty
    while Length( todo ) > 0 and c > 1 do
        j := Position(val, Minimum(val));
        g  := Remove(todo, j);
        gg := Remove(tododo, j);
        d  := Depth( g );
        f := [];

        # shift g into ind
        changed := [];
        while d < c do

            h  := ind[d];
            hh := indd[d];
            r := FactorOrder(g);
            a := LeadingExponent(g);

            # shift in
            if IsBool( h ) then
                Print("NEW\ng: ",g,"\n  gg: ",gg,"\n");
                ind[d]  := NormedPcpElement(g);
                indd[d] := gg;
                Add(f,d);
                h  := ind[d];
                hh := indd[d];
            elif not IsPrime(r) then
                b := LeadingExponent(h);
                e := Gcdex(a, b);
                if e.coeff1 <> 0 then 
                    ind[d]  := NormedPcpElement((g^e.coeff1)*(h^e.coeff2));
                    indd[d] := (gg^e.coeff1)*(hh^e.coeff2);
                    Add(f,d);
                fi;
            fi;

            # divide off
            if g = h then 
                g  := g^0;
                gg := gg^0; # Or only if gg = hh?
            else
                b  := LeadingExponent(h);
                e  := Gcdex(a,b);
                g  := g^e.coeff3 * h^e.coeff4;
                gg := gg^e.coeff3 * hh^e.coeff4;
            fi;
            d := Depth(g);
        od;

        # adjust
        c := TailLimit(ind, c);
        ReduceExpo(ind,  todo,   rels);
        #ReduceExpo(indd, tododo, rels);

        # add powers and commutators
        for d in f do
            g :=  ind[d];
            gg := indd[d];
            if rels[d] > 0 then
                r := RelativeOrderPcp(g);
                k := g ^ r;
                if Depth(k) < c then
                    Add(todo,   k);
                    Add(tododo, gg^r);
                fi;
            fi;
            for j in [1..n] do
                if not IsBool(ind[j]) then
                    k := Comm(g, ind[j]);
                    if Depth(k) < c then
                        Add(todo, k);
                        Add(tododo, Comm(gg, indd[j]));
                    fi;
                    if rels[j] = 0 then
                        k := Comm(g, ind[j]^-1);
                        if Depth(k) < c then
                            Add(todo, k);
                            Add(tododo, Comm(gg, indd[j]^-1));
                        fi;
                    fi;
                fi;
            od;
        od;

        # reduce
        sub    := Filtered([1..Length(todo)], i -> Depth(todo[i])<c);
        todo   := todo{sub};
        tododo := tododo{sub};
        val    := List(todo, x -> IGSValFun(x));
        Info(InfoPcpGrp, 3, Length(val)," versus ", ind);
    od;

    # return resulting list
    return [Filtered( ind, x -> not IsBool( x ) ),
            Filtered( indd, x -> not IsBool( x ) ) ];

    ind  := Filtered(ind,  x -> not IsBool(x));
    indd := Filtered(indd, x -> not IsBool(x));
    # CHECK_IGS here or not?
    return [ind, indd];
end );

#############################################################################
##
## IgsParallel( <gens>, <pre> )
##
InstallGlobalFunction( IgsParallel, function( gens, pre )
    return AddToIgsParallel( [], gens, [], pre );
end );

#############################################################################
##
## CgsParallel( <gens>, <pre> )
##
## parallel version of Cgs. Note: this function performes an
## induced pcs computation as well.
##
InstallGlobalFunction( CgsParallel, function( gens, pre )
    local   can,  cann,  i,  f,  e,  j,  l,  d,  r, s;

    if Length( gens ) = 0 then return []; fi;

    can  := IgsParallel( gens, pre );
    cann := can[2];
    can  := can[1];

    # first norm leading coefficients
    for i in [1..Length(can)] do
        f := NormingExponent( can[i] );
        can[i]  := can[i]^f;
        cann[i] := cann[i]^f;
    od;

    # reduce entries in matrix
    for i in [1..Length(can)] do
        e := LeadingExponent( can[i] );
        r := Depth( can[i] );
        for j in [1..i-1] do
            l := Exponents( can[j] )[r];
            if l > 0 then
                d := QuoInt( l, e );
                can[j]  := can[j]  * can[i]^-d;
                cann[j] := cann[j] * cann[i]^-d;
            elif l < 0 then
                d := QuoInt( -l, e );
                s := RemInt( -l, e );
                if s = 0 then
                    can[j] := can[j] * can[i]^d;
                    cann[j] := cann[j] * cann[i]^d;
                else
                    can[j] := can[j] * can[i]^(d+1);
                    cann[j] := cann[j] * cann[i]^(d+1);
                fi;

            fi;
        od;
    od;
    return[ can, cann ];
end );
