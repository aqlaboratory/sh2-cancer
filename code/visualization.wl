(* ::Package:: *)

(* ::Text:: *)
(*Amino Acid Color Scheme (http://www.bioinformatics.nl/~berndb/aacolour.html)*)


ColorSchemeMAEditorAAs={Black,RGBColor[7/15,13/15,8/15],RGBColor[3/5,14/15,2/5],RGBColor[1/3,11/15,1/5],RGBColor[1/3,11/15,1/5],RGBColor[3/5,3/5,1],RGBColor[7/15,13/15,8/15],RGBColor[1/3,1/3,1],RGBColor[2/5,11/15,1],RGBColor[1,4/5,7/15],RGBColor[2/5,11/15,1],RGBColor[2/5,11/15,1],RGBColor[1/3,11/15,1/5],RGBColor[14/15,2/3,2/3],RGBColor[1/3,11/15,1/5],RGBColor[1,4/5,7/15],RGBColor[1,4/15,1/3],RGBColor[1,4/15,1/3],RGBColor[2/5,11/15,1],RGBColor[3/5,3/5,1],RGBColor[3/5,3/5,1]};


(* ::Text:: *)
(*Visualizes a Position Energy Matrix (PEM). Accepts PEM object, total image size, and letter size (controls fidelty of generated image). There are a number of options. HeightFactor is a linear multiplier for all the heights. DividerThickness determines how thick the dividing line is. DividerColor is similarly self-explanatory. HighestOnTop determines whether the largest energy contributor, negative or positive, is closest to the line or furthest. So it's also "MostNegativeOnBottom". FoldOver determines whether, instead of showing the raw energies (False), the function instead shows the fold-over contribution of each letter (True). Opacity determines whether the weaker contributors (absolute value of energy) should be faded out.*)


Options[PEMPlot]={HeightFactor->1,DividerThickness->10,DividerColor->Black,HighestOnTop->False,FoldOver->False,Opacity->True,OpacityShift->0};
Options[GenerateRasterizedLetters]={ColorScheme->Rest[ColorSchemeMAEditorAAs],BackgroundFunction->(White&),TextStyle->Sequence[],Alphabet->DeleteCases[AlphaAAs,"-"]};
PEMPlot[pem_,sizeImage_/;MemberQ[{Integer,Real,Symbol},Head[sizeImage]],sizeLetter_Integer:300,opts:OptionsPattern[{PEMPlot,GenerateRasterizedLetters}]]:=Module[
	{archetypes=GenerateRasterizedLetters[sizeLetter,Sequence@@FilterRules[{opts},Options[GenerateRasterizedLetters]]],whiter,divider,orderings,blocks,blocksOrdered,heights,cols,heightsMax,heightMax},
	whiter=SetAlphaChannel[ImageResize[Graphics[],{ImageDimensions[archetypes[[1]]][[1]],1}],1];
	divider=SetAlphaChannel[ImageResize[Graphics[{},Background->OptionValue[DividerColor]],{ImageDimensions[archetypes[[1]]][[1]],OptionValue[DividerThickness]}],1];
	orderings=Transpose[Function[entries,Table[Pick[Ordering[entries],Sign[Sort[entries]],sign],{sign,{1,-1}}]]/@pem];(* Modify "orderings" to alter stacking behavior *)
	If[OptionValue[HighestOnTop],orderings=Map[Reverse,orderings,{2}]];
	heights=If[OptionValue[FoldOver],E^Abs[pem],Abs[pem]]*OptionValue[HeightFactor];(* Modify "heights" to alter heights *)
	heightMax=Max[heights];
	blocks=MapThread[
		SetAlphaChannel[ImageResize[#1,{Scaled[1],Max[#2*ImageDimensions[#1][[2]],1]}],If[OptionValue[Opacity],Min[(#2+OptionValue[OpacityShift])/heightMax,1],1]]&,
		{archetypes,#}]&/@heights;
	blocksOrdered=Replace[MapThread[{#1[[#2]]}\[Transpose]&,{blocks,#}]&/@orderings,{}->{{whiter}},{2}];
	cols=Map[ImageAssemble,blocksOrdered,{2}];
	heightsMax=Max[(ImageDimensions/@#)[[All,2]]]&/@cols;
	cols=MapThread[
		Function[{colsSection,heightMax,f},ImagePad[#,{{0,0},f[{0,heightMax-ImageDimensions[#][[2]]}]},White]&/@colsSection],
		{cols,heightsMax,{Sequence,Reverse}}];
	ImageResize[ImageAssemble[Insert[cols,ConstantArray[divider,Length[pem]],2]],sizeImage]]
PEMPlot[pem_,opts:OptionsPattern[{PEMPlot,GenerateRasterizedLetters}]]:=PEMPlot[pem,Automatic,opts]


(* ::Text:: *)
(*Takes an n*m matrix, where the entries represent the number of interactions between residues, and generates a visual representation of it.*)


Options[InteractionGraph]={EdgeWeightDownFactor->2000,EdgeColor->Black,Shape1->"Square",Shape2->"Circle",Color1->ColorSchemeBright2[[1]],Color2->ColorSchemeBright2[[4]],VertexTextStyle->{14,White,Bold,FontFamily->"Arial"},VertexSize->0.5,ImageSize->630};
InteractionGraph[mat_,OptionsPattern[]]:=Module[
	{raw=Most[ArrayRules[mat]],edges,vertexes},
	edges=UndirectedEdge[#[[1]],ToString[#[[2]]]]&/@raw[[All,1]];
	vertexes=VertexList[Graph[edges]];
	Graph[
		edges,
		EdgeStyle->Thread[Rule[edges,Directive[Thickness[#],OptionValue[EdgeColor]]&/@(raw[[All,2]]/OptionValue[EdgeWeightDownFactor])]],
		VertexCoordinates->MapIndexed[{#2[[1]],Head[#1]/.{Integer->1,String->-1}}&,vertexes],
		VertexShapeFunction->Thread[vertexes->(Head/@vertexes/.{Integer->OptionValue[Shape1],String->OptionValue[Shape2]})],
		VertexStyle->Thread[vertexes->(Head/@vertexes/.{Integer->OptionValue[Color1],String->OptionValue[Color2]})],
		VertexLabels->(#->Placed[Text[Style[#,Sequence@@OptionValue[VertexTextStyle]]],Center]&/@vertexes),
		VertexSize->OptionValue[VertexSize],ImageSize->OptionValue[ImageSize]]]
