(* ::Package:: *)

(*load the package*)
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
Get[FileNameJoin[{current, "..", "..", "AMFlow.m"}]];


(*set ibp reducer, could be "Blade", "FiniteFlow+LiteRed", "FIRE+LiteRed" or "Kira"*)
SetReductionOptions["IBPReducer" -> "FiniteFlow+LiteRed"];


(*configuration of the integral family*)
AMFlowInfo["Family"] = sunrise;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p^2 -> s};
AMFlowInfo["Propagator"] = {l1^2-msq, l2^2, (l1+l2+p)^2, (l1-p)^2, (l2-p)^2};
AMFlowInfo["Numeric"] = {msq -> 1};
AMFlowInfo["NThread"] = 4;


(*reduce target integrals to preferred master integrals*)
(*if there is no preferred master integrals, simply input with {}*)
target = {j[sunrise, 1,1,1,-3,0]};
preferred = {j[sunrise,2,1,1,0,0]};
{masters, rules} = BlackBoxReduce[target, preferred];
reduction = Thread[Keys[rules] -> Values[rules] . masters];
Put[reduction, FileNameJoin[{current, "redtable"}]];


(*construct differential equations for master integrals*)
{mastersde, vars, diffeq} = BlackBoxDiffeq[masters, {s}];
Put[{mastersde, diffeq[[1]]}, FileNameJoin[{current, "diffeq"}]];


(*perform auxiliary mass flow to obtain boundary conditions*)
AMFlowInfo["Numeric"] = {s -> 1/2, msq -> 1};
epslist = {10^-4};
sol1 = BlackBoxAMFlow[masters, epslist];
Put[sol1, FileNameJoin[{current, "sol1"}]];


(*compute at another point which can be used to perform self-consistency check*)
AMFlowInfo["Numeric"] = {s -> 1/10, msq -> 1};
epslist = {10^-4};
sol2 = BlackBoxAMFlow[masters, epslist];
Put[sol2, FileNameJoin[{current, "sol2"}]];


Quit[];
