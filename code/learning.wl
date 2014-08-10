(* ::Package:: *)

(* ::Text:: *)
(*Function takes a domain string and a peptide string and returns an input vector corresponding to all 1st and 2nd order interactions. Positions are indexed such that pY is included (even if pY terms are excluded.) FirstOrderDomain and FirstOrderPeptide both take All, None, or some list of numbers (e.g. {8,10,15}) and accordingly only return features for those positions. For 1st order terms, the pY option is irrelevant. SecondOrder takes All, None, or a list of positions, of the form {{n,m},\[Ellipsis]}, where n is the domain position and m is the peptide position. The positions asked for are the ones whose features are returned. If a position includes the pY residue, i.e. {n,10}, then depending on the pY option, a feature is either returned or not. IncludedResiduesDomain/IncludedResiduesPeptide explicitly specify the alphabet for domains and peptides. One can also choose to exclude certain residues (i.e. get set to all 0s), and that is specified by ExcludedResiduesDomain and ExcludedResiduesPeptide. The Included* options set the bandwidth of the one-hot representation.Function also accepts a list of domains and returns sum of interactions (with peptide) for those domains.*)


Options[SequenceToPositionBasedFeatures]={FirstOrderDomain->All,FirstOrderPeptide->All,SecondOrder->All,pY->False,Structured->False,IncludedResiduesDomain->AlphaAAs,IncludedResiduesPeptide->AlphaAAs,ExcludedResiduesDomain->{},ExcludedResiduesPeptide->{},PYResidue->"y"};
SequenceToPositionBasedFeatures[strDomain_String,strPeptide_String,OptionsPattern[]]:=Module[
	{features1stOrderDomain,features1stOrderPeptide,features2ndOrder,charsDomain,charsPeptide,rules1stOrderDomain,rules1stOrderPeptide,rules2ndOrder,
	opts1stOrderDomain=OptionValue[FirstOrderDomain]/.None->{},opts1stOrderPeptide=OptionValue[FirstOrderPeptide]/.None->{},opts2ndOrder=OptionValue[SecondOrder],optsPY=OptionValue[pY],optsStructured=OptionValue[Structured],pY=OptionValue[PYResidue],alphaIncludedDomain=OptionValue[IncludedResiduesDomain],alphaIncludedPeptide=OptionValue[IncludedResiduesPeptide],alphaExcludedDomain=OptionValue[ExcludedResiduesDomain],alphaExcludedPeptide=OptionValue[ExcludedResiduesPeptide]},
	rules1stOrderDomain=Dispatch[Join[
		Thread[alphaIncludedDomain->IdentityMatrix[Length[alphaIncludedDomain]]],
		Thread[alphaExcludedDomain->ConstantArray[0,{Length[alphaExcludedDomain],Length[alphaIncludedDomain]}]]]];
	rules1stOrderPeptide=Dispatch[Join[
		Thread[alphaIncludedPeptide->IdentityMatrix[Length[alphaIncludedPeptide]]],
		{pY->{}},
		Thread[alphaExcludedPeptide->ConstantArray[0,{Length[alphaExcludedPeptide],Length[alphaIncludedPeptide]}]]]];
	rules2ndOrder=With[
		{tups=Flatten[Outer[Tuples[{##}]&,{alphaIncludedDomain,alphaExcludedDomain},{alphaIncludedPeptide,alphaExcludedPeptide},1],1]},
		Dispatch[Join[
			Thread[First[tups]->IdentityMatrix[Length[First[tups]]]],
			Flatten[Thread[#->ConstantArray[0,{Length[#],Length[First[tups]]}]]&/@Rest[tups],1],
			If[optsPY,
				Join[
					Thread[Thread[{alphaIncludedDomain,pY}]->IdentityMatrix[Length[alphaIncludedDomain]]],
					Thread[Thread[{alphaExcludedDomain,pY}]->ConstantArray[0,{Length[alphaExcludedDomain],Length[alphaIncludedDomain]}]]],
				{{_,pY}->{}}]]]];
	{charsDomain,charsPeptide}=Characters/@{strDomain,strPeptide};
	features1stOrderDomain=If[optsStructured,Sequence,Flatten][charsDomain[[opts1stOrderDomain]]/.rules1stOrderDomain];
	features1stOrderPeptide=If[optsStructured,Sequence,Flatten][charsPeptide[[opts1stOrderPeptide]]/.rules1stOrderPeptide];
	features2ndOrder=Outer[List,charsDomain,charsPeptide];
	features2ndOrder=If[optsStructured,Sequence,Flatten][Switch[opts2ndOrder,
		All,Flatten[features2ndOrder,1],
		None,{},
		_,Extract[features2ndOrder,opts2ndOrder]]/.rules2ndOrder];
	If[optsStructured,List,Join][features1stOrderDomain,features1stOrderPeptide,features2ndOrder]]
SequenceToPositionBasedFeatures[strDomain_List,strPeptide_String,opts:OptionsPattern[]]:=Total[SequenceToPositionBasedFeatures[#,strPeptide,opts]&/@strDomain]


(* ::Text:: *)
(*Takes a domain string and a peptide string as well as sets of indices of interacting pairs, and generates an input vector corresponding to all the (2nd order) interactions present in the indices given. The way this function works is residue-based, i.e. each feature corresponds to the count of the number of times an interacting residue pair is observed, regardless of position.*)


SequenceToResidueBasedFeatures[strDomain_String,strPeptide_String,posGroups_]:=Module[
	{features,charsDomain,charsPeptide,tuples,counts,
	alphaAAs=AlphaAAs},
	tuples=Join[Tuples[alphaAAs,2],({#,"y"}&/@alphaAAs)];
	{charsDomain,charsPeptide}=Characters/@{strDomain,strPeptide};
	features=Outer[List,charsDomain,charsPeptide];
	counts=Switch[Depth[posGroups],
		4,Append[Thread[Rule@@Transpose[Tally[Extract[features,#]]]],{_,_}->0]&/@posGroups,
		3,{{_,_}->0}];
	Flatten[(tuples/.#)&/@counts]]


(* ::Text:: *)
(*Takes a contact map and a cutoff distance, and returns a list of interactions that are shorter than the cutoff distance and that are grouped by adjacency in terms of their positions (two interactions are adjacent if both the domain and peptide residues (in the two interactions) are adjacent physically.)*)


ContactMapToConnectedPositions[matContacts_,cutoff_]:=Module[{pos,amat},
	pos=Position[matContacts,d_/;d<cutoff];
	amat=Outer[ChessboardDistance,pos,pos,1]/.(d_/;d>1)->0;
	pos[[#]]&/@ConnectedComponents[AdjacencyGraph[amat]]]
