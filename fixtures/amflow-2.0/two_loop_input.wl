(* ::Package:: *)

(*load the package*)
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
Get[FileNameJoin[{current, "..", "..", "AMFlow.m"}]];


(*set ibp reducer, could be "Blade", "FiniteFlow+LiteRed", "FIRE+LiteRed" or "Kira"*)
SetReductionOptions["IBPReducer" -> "FiniteFlow+LiteRed"];


(*read results from cached files whenever possible, False by default*)
SetAMFOptions["UseCache" -> True];


(*configuration of the integral family*)
AMFlowInfo["Family"] = tt;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1-p2-p3};
AMFlowInfo["Replacement"] = {p1^2 -> 0, p2^2 -> 0, p3^2 -> msq, p4^2 -> msq, (p1+p2)^2 -> s, (p1+p3)^2 -> t};
AMFlowInfo["Propagator"] = {l1^2, (l1+p1)^2, (l1+p1+p2)^2, l2^2, -msq+(l2+p3)^2, (l2+p3+p4)^2, (l1+l2)^2, (l1-p3)^2, (l2+p1)^2};
AMFlowInfo["Numeric"] = {s -> 30, t -> -10/3, msq -> 1};
AMFlowInfo["NThread"] = 4;


(*SolveIntegrals: computes given integrals with given precision goal up to given eps order*)
(*returned is a list of replacement rules like {j1 \[Rule] v1, j2 \[Rule] v2, ...}, where j1, j2, ... are integrals and v1, v2, ... are their results*)
target = {j[tt,1,1,1,1,1,1,1,-3,0], j[tt,1,1,1,1,1,1,1,-2,-1], j[tt,1,1,1,1,1,1,1,-1,-2], j[tt,1,1,1,1,1,1,1,0,-3]};
precision = 20;
epsorder = 4;
auto = SolveIntegrals[target, precision, epsorder];


(*save the results*)
Put[auto, FileNameJoin[{current, "auto"}]];


(*manual computation*)
(*suggested configuration in SolveIntegrals, can be regenerated using GenerateNumericalConfig*)
(*of course users can define their own configurations without using this function*)
{epslist, workingpre, xorder, exorder, leading} = GenerateNumericalConfig[precision, epsorder];


(*set options for auxiliary mass flow*)
(*"WorkingPre": working precision for numerical evaluations, 100 by default*)
(*"XOrder": truncated order of eta expansion during the evaluation, 100 by default*)
(*"ExtraXOrder": extra order for some internal algorithms, 30 by default*)
SetAMFOptions["WorkingPre" -> workingpre, "XOrder" -> xorder, "ExtraXOrder" -> exorder];


(*compute target integrals at give eps values*)
(*structure of soleps: {int1 \[Rule] {int1[eps1], int1[eps2], ...}, int2 \[Rule] {int2[eps1], int2[eps2], ...}, ...}*)
soleps = BlackBoxAMFlow[target, epslist];


(*fit the expansion of target integrals*)
exp = FitEps[epslist, #, leading]&/@Values[soleps];
man = Thread[Keys[soleps] -> Chop@N[Normal@Series[exp, {eps, 0, leading+epsorder}], precision]];


(*save the results*)
Put[man, FileNameJoin[{current, "man"}]];


Quit[];
