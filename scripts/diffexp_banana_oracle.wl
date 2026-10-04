(* Independently written driver for the pinned equal-mass banana example.
   Required environment: DIFFEXP_ROOT (pinned checkout), DIFFEXP_OUTPUT
   (a NEW output directory). The upstream checkout is loaded unchanged.
   Use one process/thread and a bound, for example:
     OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     timeout --signal=TERM --kill-after=10s 600s \
     bern-wolfram -noinit -noprompt -script scripts/diffexp_banana_oracle.wl
   All four masters through epsilon^4 are initialized from the notebook's
   analytic partial gamma boundary. No numerical reference seeds are used.
   Baseline settings match the notebook; the second run increases precision
   and order. Saved Pade values at t=20 are retained as additional raw output.
*)
$HistoryLength=0;
driverFile=ExpandFileName[$InputFileName];
fail[m_]:=(Print["ORACLE FAILURE: ",m];Exit[1]);
requiredEnvironment[name_]:=Module[{value=Environment[name]},
 If[!StringQ[value]||StringLength[value]==0,fail["Set "<>name<>" before running"]];
 ExpandFileName[value]];
base=requiredEnvironment["DIFFEXP_OUTPUT"];
root=requiredEnvironment["DIFFEXP_ROOT"];
str[v_]:=ToString[v,InputForm];sha[f_]:=IntegerString[FileHash[f,"SHA256"],16,64];
scalar[v_]:=<|"real_wolfram"->str[Re[v]],"imaginary_wolfram"->str[Im[v]],"precision_wolfram"->str[Precision[v]],"accuracy_wolfram"->str[Accuracy[v]]|>;
If[sha[FileNameJoin[{root,"DiffExp.m"}]]=!="67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7",fail["Pinned package mismatch"]];
If[FileExistsQ[base]||DirectoryQ[base],fail["DIFFEXP_OUTPUT must be a new directory"]];
If[CreateDirectory[base,CreateIntermediateDirectories->True]===$Failed,fail["Cannot create output directory"]];
matrixHashes={"d87a9039ab70226a0afbb4c4c2c37f2fbfde6c638c30a2771087734ae982ef1b","c4c597adeda5872d532d3edb82a11bdfbd04456f39e6bd1d5fd538390f96fb71","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab"};
Do[file=FileNameJoin[{root,"Examples","Banana_EqualMass_Matrices","dt_"<>ToString[k]<>".m"}];
 If[!FileExistsQ[file]||sha[file]=!=matrixHashes[[k+1]],fail["Pinned matrix missing or hash mismatch"]],{k,0,4}];
Get[FileNameJoin[{root,"DiffExp.m"}]];
rawBoundary={"?","?",\[Epsilon]*(1+3*\[Epsilon])*(1+4*\[Epsilon])*E^(3*EulerGamma*\[Epsilon])*(-4*Gamma[\[Epsilon]]^3/t+6*(-1/t)^(1+\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^2*Gamma[\[Epsilon]]^3/Gamma[-2*\[Epsilon]]+8*(-1/t)^(1+2*\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^3*Gamma[\[Epsilon]]*Gamma[2*\[Epsilon]]/Gamma[-3*\[Epsilon]]+3*(-1/t)^(1+3*\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^4*Gamma[3*\[Epsilon]]/Gamma[-4*\[Epsilon]]),E^(3*EulerGamma*\[Epsilon])*\[Epsilon]^3*Gamma[\[Epsilon]]^3};
records={};
CheckAbort[
 Do[
  wp=setting[[1]];order=setting[[2]];
  DiffExp`LoadConfiguration[{DiffExp`MatrixDirectory->FileNameJoin[{root,"Examples","Banana_EqualMass_Matrices"}],DiffExp`EpsilonOrder->4,WorkingPrecision->wp,DiffExp`ChopPrecision->Floor[wp/2],DiffExp`ExpansionOrder->order,DiffExp`DivisionOrder->3,DiffExp`UseMobius->True,DiffExp`UsePade->True,DiffExp`DeltaPrescriptions->{t-16+I*\[Delta]},DiffExp`Verbosity->1,"Parallel"->False}];
  preparation=AbsoluteTiming[DiffExp`PrepareBoundaryConditions[rawBoundary,{t->-1/x}]];
  Print["STAGE: infinity initialization WP",wp," order",order];
  initial=AbsoluteTiming[DiffExp`TransportTo[preparation[[2]],{t->-1}]];
  If[Dimensions[initial[[2,2]]]=!={4,5}||!AllTrue[Flatten[initial[[2,2]]],NumberQ],fail["Malformed t=-1 result"]];
  Put[initial[[2]],FileNameJoin[{base,"minus1-order"<>ToString[order]<>".wl"}]];
  Print["STAGE: physical continuation WP",wp," order",order];
  physical=AbsoluteTiming[DiffExp`TransportTo[initial[[2]],{t->x},32,True]];
  endpoint=physical[[2,1]];
  If[Dimensions[endpoint[[2]]]=!={4,5}||!AllTrue[Flatten[endpoint[[2]]],NumberQ],fail["Malformed t=32 result"]];
  Put[physical[[2]],FileNameJoin[{base,"physical-order"<>ToString[order]<>".wl"}]];
  saved=AbsoluteTiming[DiffExp`ToPiecewise[physical[[2]],True]];
  at20=Table[saved[[2,i,j]][20],{i,4},{j,5}];
  If[!AllTrue[Flatten[at20],NumberQ],fail["Malformed saved-series t=20 result"]];
  AppendTo[records,<|"working_precision"->wp,"chop_precision"->Floor[wp/2],"expansion_order"->order,"preparation_seconds"->preparation[[1]],"infinity_to_minus1_seconds"->initial[[1]],"minus1_to_32_seconds"->physical[[1]],"saved_pade_construction_seconds"->saved[[1]],"minus1_values"->Map[scalar,initial[[2,2]],{2}],"endpoint32_values"->Map[scalar,endpoint[[2]],{2}],"saved20_values"->Map[scalar,at20,{2}],"minus1_errors"->Map[str,initial[[2,3]],{2}],"endpoint32_errors"->Map[str,endpoint[[3]],{2}]|>];
  Export[FileNameJoin[{base,"runs-so-far.json"}],records,"RawJSON"],
 {setting,{{500,50},{600,75}}}];
 errors=Table[
   lo=If[kind==="minus1",Get[FileNameJoin[{base,kind<>"-order50.wl"}]][[2]],Get[FileNameJoin[{base,"physical-order50.wl"}]][[1,2]]];
   hi=If[kind==="minus1",Get[FileNameJoin[{base,kind<>"-order75.wl"}]][[2]],Get[FileNameJoin[{base,"physical-order75.wl"}]][[1,2]]];
   Max[Abs[Flatten[lo-hi]]],{kind,{"minus1","physical"}}];
 If[!TrueQ[Max[errors]<10^-20],fail[{"Precision/order refinement below 20 digits",errors}]];
 exported=Export[FileNameJoin[{base,"result.json"}],<|"status"->"passed","upstream_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd","package_sha256"->sha[FileNameJoin[{root,"DiffExp.m"}]],"driver_sha256"->sha[driverFile],"wolfram_version"->$Version,"settings"-><|"epsilon_order"->4,"dimension"->"2-2*epsilon","division_order"->3,"mobius"->True,"pade"->True,"parallel"->False,"prescription"->"t-16+I*delta"|>,"runs"->records,"max_absolute_refinement_errors_minus1_32"->(str/@errors),"asserted_digits"->20,"scope"->"Complete equal-mass banana four-master DE, epsilon orders0..4, analytic partial infinity boundary. No IBP reduction timing.","notebook_deviations"->{"Verbosity1 instead of2; second500→600 digit/order50→75 refinement added."}|>,"RawJSON"];
 If[exported===$Failed,fail["Failed result export"]];
 Print["ORACLE PASS ",str[errors]];Exit[0],fail["DiffExp aborted"]];
