(* ::Package:: *)

(* ::Text:: *)
(*Computes recall, precision, Subscript[F, 1], and accuracy between two vectors of 0s and 1s.*)


RPFA[exp_,pred_]:=Module[
	{npred=1-pred,nexp=1-exp,tp,fp,fn,recall,precision,f1,accuracy},
	tp=pred.exp;
	fp=pred.nexp;
	fn=npred.exp;
	recall=Quiet[tp/(tp+fn),{Power::"infy",Infinity::"indet"}]/.Indeterminate->0;
	precision=Quiet[tp/(tp+fp),{Power::"infy",Infinity::"indet"}]/.Indeterminate->0;
	f1=Quiet[(2 precision recall)/(precision+recall),{Power::"infy",Infinity::"indet"}]/.Indeterminate->0;
	accuracy=1-HammingDistance[exp,pred]/Length[exp];
	{recall,precision,f1,accuracy}]


(* ::Text:: *)
(*Computes true positive rate and false positive rate. Same inputs as RPFA.*)


TFPR[exp_,pred_]:=Module[{npred=1-pred,nexp=1-exp,tp,fp,tn,fn,tpr,fpr},
	tp=pred.exp;
	fp=pred.nexp;
	tn=npred.nexp;
	fn=npred.exp;
	tpr=tp/(tp+fn);
	fpr=fp/(fp+tn);
	{fpr,tpr}]


(* ::Text:: *)
(*Generates an ROC curve. If the number of steps requested is smaller than the number of distinct entries, then the returned list is of length equal to the number of distinct entries. If \[Infinity] is used for nSteps, the function uses all distinct entries.*)


ROC[exp_,pred_,nSteps_:100]:=Module[
	{steps=DeleteDuplicates[Sort[pred]],len},
	len=Length[steps];
	steps=steps[[Append[Range[1,len,Max[\[LeftCeiling]len/(nSteps-1)\[RightCeiling],1]],len]]];
	Table[TFPR@@Unitize[Chop[{exp,pred},cutoff]],{cutoff,steps}]]


(* ::Text:: *)
(*Returns Area Under the Curve from an ROC curve.*)


AUC[roc_]:=Integrate[Interpolation[DeleteDuplicates[roc,(#1[[1]]===#2[[1]])&],InterpolationOrder->1][x],{x,0,1}]
