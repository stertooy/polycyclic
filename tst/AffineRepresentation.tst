gap> START_TEST( "AffineRepresentation.tst" );

# Terminate at the trivial group and after the final polycyclic generator.
gap> G := AbelianPcpGroup( [] );;
gap> AffineRepresentation( G );
[  ]
gap> G := AbelianPcpGroup( 1 );;
gap> mats := AffineRepresentation( G );;
gap> IsMatrixRepresentation( G, GeneratorsOfPcp( Pcp( G ) ), mats );
true
gap> GroupHomomorphismByImages( G, Group( mats ),
>        GeneratorsOfPcp( Pcp( G ) ), mats ) <> fail;
true

# Intermediate quotients can have more group generators than Pcp generators.
gap> G := UnitriangularPcpGroup( 3, 0 );;
gap> mats := AffineRepresentation( G );;
starting level 2
  got cohomology with orders [ 0, 0 ]
starting level 3
  got cohomology with orders [ 0, 0, 0 ]
gap> IsMatrixRepresentation( G, GeneratorsOfPcp( Pcp( G ) ), mats );
true
gap> GroupHomomorphismByImages( G, Group( mats ),
>        GeneratorsOfPcp( Pcp( G ) ), mats ) <> fail;
true

# A subgroup can have fewer Pcp generators than its ambient collector.
gap> H := Subgroup( G, [ G.2, G.3 ] );;
gap> mats := AffineRepresentation( H );;
starting level 2
  got cohomology with orders [ 0, 0 ]
gap> IsMatrixRepresentation( H, GeneratorsOfPcp( Pcp( H ) ), mats );
true
gap> GroupHomomorphismByImages( H, Group( mats ),
>        GeneratorsOfPcp( Pcp( H ) ), mats ) <> fail;
true
gap> H := Subgroup( G, [ G.1^2, G.2, G.3 ] );;
gap> mats := AffineRepresentation( H );;
starting level 2
  got cohomology with orders [ 0, 0 ]
starting level 3
  got cohomology with orders [ 0, 0, 0 ]
gap> IsMatrixRepresentation( H, GeneratorsOfPcp( Pcp( H ) ), mats );
true
gap> GroupHomomorphismByImages( H, Group( mats ),
>        GeneratorsOfPcp( Pcp( H ) ), mats ) <> fail;
true

gap> STOP_TEST( "AffineRepresentation.tst" );
