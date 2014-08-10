(* ::Package:: *)

(* ::Text:: *)
(*Generates a Position Energy Matrix (PEM) object from an SH2 domain sequence and the \[Beta] coefficients as outputted by ImportBetaMatrix (single \[Lambda]). Options exist to turn off the first-order domain and peptide energies if desired.*)


Options[PEM]={YPosition->7,AlphabetDomain->DeleteCases[AlphaAAs,"-"],FirstOrderDomainEnergies->False,FirstOrderPeptideEnergies->True};
PEM[seqSH2_,\[Beta]sRaw_,OptionsPattern[]]:=Block[
	{\[Beta]s=MapAt[Transpose[#,{1,3,2,4}]&,\[Beta]sRaw,3],
	 alphaAADomain=OptionValue[AlphabetDomain],posY=OptionValue[YPosition],dummy,rules},
	dummy=SparseArray[{},Dimensions[\[Beta]s[[2]]]];
	rules=Dispatch[Append[Thread[alphaAADomain->Range[Length[alphaAADomain]]],"-"->0]];
	ReplacePart[
		Normal[Total[MapThread[If[#3=!=0,#2[[#3]]+If[OptionValue[FirstOrderDomainEnergies],#1[[#3]],0],dummy]&,{\[Beta]s[[1]],\[Beta]s[[3]],Characters[seqSH2]/.rules}]]],
		{posY,-1}->0]+If[OptionValue[FirstOrderPeptideEnergies],\[Beta]s[[2]],0]]


(* ::Text:: *)
(*The functions below compute, given as input the metadata (with indexes to the matrix of the individual SH2-phosphosites probabilities) of SH2 domains and peptides, the direct (protein-level) interaction probabilities and the perturbation probabilities, respectively. The perturbation functions return a list of probabilities, two to be precise, the first for going from inactive to active, and the second from active to inactive. The total perturbation probability is just the sum of the two since these two outcomes are mutually exclusive. In the case that ModelPosterior is set to False, the list is just one element, since the priors do not distinguish between those two outcomes.*)
(**)
(*SystemPosterior controls whether the protein context is taken into account. ModelPosterior controls where the domain model probabilities are used or not. Model controls the semantics between the "V1" and "V2" models. "V2" is what is described in AlQuraishi et. al 2014. What is returned is slightly different depending on the model. For "V1", the states corresponding to the probabilities of each single edge being active are returned, but nothing else, because the edges are independent and so one can compute anything from that. For "V2", the function returns the entire ensemble, i.e. every allowable state with its probability. These "V2" probabilities are all unnormalized, i.e. they need to be divided by the partition function (one can say the same about "V1" but that's vacuous).*)
(**)
(*This function supports multiple experiments/predictions per entry (peptide/domain, and so a gene containing multiple peptides and multiple domains will further decompose into multiple entries for each peptide/domain combo in the gene). How to choose among these entries is specificed by the PoolingFunction option.*)


DistributeDefinitions[DirectedPerturbationProbability,DirectedInteractionProbability,DirectedPerturbationProbabilities,DirectedInteractionProbabilities,DirectedInteractionProbabilityDecomposition];


Options[DirectedInteractionProbabilities]={Model->"V1",PoolingFunction->Max};
DirectedInteractionProbabilities[{rowsProts_,colsProts_},matProbs_,opts:OptionsPattern[]]:=Outer[DirectedInteractionProbability[{##},matProbs,opts]&,rowsProts,colsProts,1]
DirectedInteractionProbability[{rowsProt_,colsProt_},matProbs_,opts:OptionsPattern[DirectedInteractionProbabilities]]:=Block[
	{probs=DirectedInteractionProbabilityDecomposition[{rowsProt,colsProt},matProbs,opts][[All,2]]},
	Switch[OptionValue[Model],"V1",1-(Times@@(1-probs)),"V2",1/(1+Last[probs]/Total[Most[probs]])]]
DirectedInteractionProbabilityDecomposition[{rowsProt_,colsProt_},matProbs_,OptionsPattern[DirectedInteractionProbabilities]]:=Block[
	{probs,probsStructured,dims,statesMats,statesProds,fPooling=OptionValue[PoolingFunction]},
	probs=fPooling[Extract[matProbs,Tuples[#]]/.{}->{0}]&/@Tuples[Map[If[Head[#]===Integer,{#},#]&,{rowsProt[[All,1]],colsProt[[All,1]]},{2}]];
	Switch[OptionValue[Model],
	"V1",
		statesMats=Partition[#,Length[colsProt]]&/@Array[UnitVector[Length[probs],#]&,Length[probs]];
		Thread[statesMats->probs],
	"V2",
		probsStructured=Partition[probs,Length[colsProt]];
		dims=Dimensions[probsStructured];
		statesMats=Permutations[IdentityMatrix[dims[[2]]]~Join~ConstantArray[0,dims],{dims[[1]]}];
		statesProds=Times@@Flatten[#[[1]]*probsStructured+#[[2]]]&/@Transpose[{(2statesMats-1),1-statesMats}];
		Thread[statesMats->statesProds]]]


Options[DirectedPerturbationProbabilities]={ModelPosterior->True,SystemPosterior->True,Parallelize->False};
DirectedPerturbationProbabilities[{rowsProts_,colsProts_},matProbs_,rowProtsPerturbedPost_,opts:OptionsPattern[{DirectedPerturbationProbabilities,DirectedInteractionProbabilities}]]:=With[
	{posRowProtsPerturbedPost=If[OptionValue[Parallelize],Parallelize,Sequence][Position[rowsProts[[All,All,2;;3]],#,{2},1][[1]]&/@rowProtsPerturbedPost[[All,2;;3]]]},
	SparseArray[If[OptionValue[Parallelize],Parallelize,Sequence][MapThread[
		Function[
			{rowProtPerturbedPost,posRowProtPerturbedPost},
			SparseArray[ReplacePart[
				ConstantArray[If[OptionValue[ModelPosterior],{0,0},{0}],Length/@{rowsProts,colsProts}],
				posRowProtPerturbedPost[[1]]->(DirectedPerturbationProbability[{rowsProts[[posRowProtPerturbedPost[[1]]]],#},matProbs,rowProtPerturbedPost,posRowProtPerturbedPost[[2]],opts]&/@colsProts)]]],
		{rowProtsPerturbedPost,posRowProtsPerturbedPost}]]]]
DirectedPerturbationProbability[{rowsProt_,colsProt_},matProbs_,rowProtPerturbedPost_,opts:OptionsPattern[{DirectedPerturbationProbabilities,DirectedInteractionProbabilities}]]:=DirectedPerturbationProbability[{rowsProt,colsProt},matProbs,rowProtPerturbedPost,Position[rowsProt[[All,2;;3]],rowProtPerturbedPost[[2;;3]]]/.{{idx_}}:>idx,opts]
DirectedPerturbationProbability[{rowsProt_,colsProt_},matProbs_,rowProtPerturbedPost_,pos_,opts:OptionsPattern[{DirectedPerturbationProbabilities,DirectedInteractionProbabilities}]]:=Module[
	{rowsProtUnperturbed,rowProtPerturbedPre,rowsProtPre,rowsProtPost,pPriorActive,pPreActive,pPostActive,pPreInactive,pPostInactive,decompPre,decompPost,stateInactive,statesActive,
	 optsDIP=FilterRules[{opts},Options[DirectedInteractionProbabilities]],
	 optModelPosterior=OptionValue[ModelPosterior],
	 optSystemPosterior=OptionValue[SystemPosterior],
	 optModel=OptionValue[Model]},
	If[pos==={},
		If[optModelPosterior,{0,0},{0}],
		If[(\[Not]optModelPosterior\[And](\[Not]optSystemPosterior\[Or]optModel==="V2")),
			{1},
			If[optModel==="V1"\[Or]\[Not]optSystemPosterior,
				rowsProtUnperturbed=Delete[rowsProt,pos];
				rowProtPerturbedPre=Extract[rowsProt,pos];
				If[optSystemPosterior,pPriorActive=DirectedInteractionProbability[{rowsProtUnperturbed,colsProt},matProbs,Sequence@@optsDIP]];
				If[optModelPosterior,
					pPreActive=DirectedInteractionProbability[{{rowProtPerturbedPre},colsProt},matProbs,Sequence@@optsDIP];
					pPostActive=DirectedInteractionProbability[{{rowProtPerturbedPost},colsProt},matProbs,Sequence@@optsDIP]];
					If[optSystemPosterior,1-pPriorActive,1]*If[optModelPosterior,{(1-pPreActive)pPostActive,pPreActive(1-pPostActive)},{1}],

				rowsProtPre=rowsProt;
				rowsProtPost=ReplacePart[rowsProt,pos->rowProtPerturbedPost];
				{decompPre,decompPost}=Table[DirectedInteractionProbabilityDecomposition[{rsP,colsProt},matProbs,Sequence@@optsDIP],{rsP,{rowsProtPre,rowsProtPost}}];
				stateInactive=ConstantArray[0,Length/@{rowsProt,colsProt}];
				statesActive=Table[ReplacePart[stateInactive,{pos,j}->1],{j,Length[colsProt]}];
				{{pPreInactive,pPostInactive},{pPreActive,pPostActive}}=Outer[Total[#1/.#2]/Total[#2[[All,2]]]&,{{stateInactive},statesActive},{decompPre,decompPost},1];
				{pPreInactive*pPostActive,pPreActive*pPostInactive}]]]]


(* ::Text:: *)
(*Undirected versions of the Interaction and Perturbation Probabiiities functions.*)


UndirectedInteractionProbabilities[{rowsProts_,colsProts_},matProbs_,opts:OptionsPattern[DirectedInteractionProbabilities]]:=Module[
	{labelsRowProts,labelsColProts,labelsMerged,mapRowProts,mapColProts,rowsProtsPopped,colsProtsPopped,prots},
	{labelsRowProts,labelsColProts}={rowsProts[[All,1,2]],colsProts[[All,1,2]]};
	labelsMerged=Union[labelsRowProts,labelsColProts];
	{mapRowProts,mapColProts}=Map[Position[labelsMerged,#][[1,1]]&,{labelsRowProts,labelsColProts},{2}];
	{rowsProtsPopped,colsProtsPopped}=MapThread[ReplacePart[ConstantArray[{},Length[labelsMerged]],Thread[#1->#2]]&,{{mapRowProts,mapColProts},{rowsProts,colsProts}}];
	prots={rowsProtsPopped,colsProtsPopped}\[Transpose];
	Outer[UndirectedInteractionProbability[{##},matProbs,opts]&,prots,prots,1]]
UndirectedInteractionProbability[{prot1_,prot2_},matProbs_,opts:OptionsPattern[DirectedInteractionProbabilities]]:=Module[
	{interactions=Cases[{{prot1[[1]],prot2[[2]]},{prot2[[1]],prot1[[2]]}},{{__},{__}}],probs},
	probs=DirectedInteractionProbability[#,matProbs,opts]&/@interactions;
	1-(Times@@(1-probs))]
