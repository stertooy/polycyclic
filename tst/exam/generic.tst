#
# UnitriangularPcpGroup
#

#
gap> G:=UnitriangularPcpGroup(4,4);
fail
gap> G:=UnitriangularPcpGroup(0,2);
fail

#
gap> G:=UnitriangularPcpGroup(4,3);
Pcp-group with orders [ 3, 3, 3, 3, 3, 3 ]
gap> IsConfluent(Collector(G));
true
gap> hom:=GroupHomomorphismByImages( G, Group(G!.mats), Igs(G), G!.mats);;
gap> IsBijective(hom);
true

#
gap> G:=UnitriangularPcpGroup(5,0);
Pcp-group with orders [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> IsConfluent(Collector(G));
true
gap> hom:=GroupHomomorphismByImages( G, Group(G!.mats), Igs(G), G!.mats);;

#
gap> fnc1:={i,j}->PcpGroupToPcGroup(BlowUpPcpPGroup(PcGroupCode(i,j)));;
gap> fnc2:=G->[Size(G),CodePcGroup(G)];;
gap> fnc:={i,j}->fnc2(fnc1(i,j));;
gap> fnc(0,2);
[ 2, 0 ]
gap> fnc(5,4);
[ 8, 291 ]
gap> fnc(0,4);
[ 8, 0 ]
gap> List([323,33,36,2343,0], i->fnc(i,8));
[ [ 128, 4648286596389404735 ], [ 128, 283691314053165 ], [ 128, 334923589489112660889148235683641060 ],
  [ 128, 89905366469666764161151872091306780267034367 ], [ 128, 0 ] ]
gap> fnc(0,3);
[ 9, 0 ]
gap> BlowUpPcpPGroup(PcGroupCode(5,9));
Pcp-group with orders [ 3, 3, 3, 3, 3, 3, 3, 3 ]
gap> BlowUpPcpPGroup(PcGroupCode(0,9));
Pcp-group with orders [ 3, 3, 3, 3, 3, 3, 3, 3 ]
