(* Execute the original GPL notebook helper definitions directly from the pinned
   checkout; this harness does not copy them into the Rust implementation. *)
$HistoryLength=0;
root=Environment["DIFFEXP_ROOT"];
If[!StringQ[root]||root==="",Print["Set DIFFEXP_ROOT to the pinned checkout"];Exit[1]];
base=Environment["MPL_OUTPUT"];
kind=Environment["MPL_CASE"];
order=ToExpression[Environment["MPL_ORDER"]];
fail[m_]:=(Print["ORACLE FAILURE: ",m];Exit[1]);
sha[p_]:=IntegerString[FileHash[p,"SHA256"],16,64];
package=FileNameJoin[{root,"DiffExp.m"}];
notebook=FileNameJoin[{root,"Examples","MultiplePolylogarithms.nb"}];
If[sha[package]=!="67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7",fail["package hash"]];
If[!MemberQ[{{"trailing",50},{"trailing",75},{"complex",50},{"complex",75},{"weight20",100}},{kind,order}],fail["profile"]];
If[!StringQ[base]||base===""||!DirectoryQ[base],fail["MPL_OUTPUT must name a fresh existing directory"]];
If[FileExistsQ[FileNameJoin[{base,"result.json"}]],fail["existing result"]];
If[sha[notebook]=!="028003cce75bdfb911d3bbf8cfc313bb6cb5454d8ce44722c4eba8bfcaca24ba",fail["notebook hash"]];
Get[package];
boxes=Cases[Get[notebook],Cell[BoxData[b_],"Input",___]:>b,Infinity];
If[Length[boxes]=!=7,fail["notebook shape"]];
Scan[Function[i,Scan[ReleaseHold,Flatten[{ToExpression[boxes[[i]],StandardForm,HoldComplete]}]]],{2,4}];
TmpDirectory=CreateDirectory[FileNameJoin[{base,"matrices"}]];
DiffExpConfiguration=<|EpsilonOrder->0,WorkingPrecision->250,ChopPrecision->225,ExpansionOrder->order,UseMobius->True,UsePade->False,DivisionOrder->3,Verbosity->1,"Parallel"->False|>;
args=Switch[kind,"trailing",{1,-10,0,4},"complex",{10,-10+I,-1/2,-50,1},"weight20",Range[21]];
scalar[v_]:=<|"real_wolfram"->ToString[Re[v],InputForm],"imaginary_wolfram"->ToString[Im[v],InputForm],"precision_wolfram"->ToString[Precision[v],InputForm],"accuracy_wolfram"->ToString[Accuracy[v],InputForm],"real_accuracy_wolfram"->ToString[Accuracy[Re[v]],InputForm],"imaginary_accuracy_wolfram"->ToString[Accuracy[Im[v]],InputForm],"real_precision_wolfram"->ToString[Precision[Re[v]],InputForm],"imaginary_precision_wolfram"->ToString[Precision[Im[v]],InputForm]|>;
CheckAbort[
 Print["PROFILE START ",kind," order ",order];
 measured=AbsoluteTiming[Apply[G,args]/.G->GEvaluate];
 value=Coefficient[measured[[2]],pm,0];
 uncertainty=Coefficient[measured[[2]],pm];
 If[!NumberQ[value]||!NumberQ[uncertainty]||!TrueQ[Precision[value]>=20],fail[{"malformed or insufficient recorded precision",measured}]];
 Put[measured,FileNameJoin[{base,"result.wl"}]];
 result=<|"status"->"passed","case"->kind,"order"->order,"arguments_wolfram"->ToString[args,InputForm],"value"->scalar[value],"reported_uncertainty_wolfram"->ToString[uncertainty,InputForm],"wall_seconds_original_wrapper"->measured[[1]],"settings"-><|"working_precision"->250,"chop_precision"->225,"expansion_order"->order,"epsilon_order"->0,"division_order"->3,"mobius"->True,"pade"->False,"parallel"->False,"verbosity"->1|>,"upstream_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd","source_sha256"->sha[package],"notebook_sha256"->sha[notebook],"driver_sha256"->sha[$InputFileName],"wolfram_version"->$Version,"scope"->"Original notebook GEvaluate workload, including its shuffle/logarithm extraction and all underlying transports. Native direct triangular-DE timing is different work; no matched algorithm ratio.","notebook_deviations"->{"Verbosity1 preserves logs; original0.","Pinned local source and target-owned temporary matrix directory; single-process mode."},"validation_note"->"This status validates a finite returned value with >=20 recorded precision; independent native/refinement comparison is a separate pending gate."|>;
 Export[FileNameJoin[{base,"result.json"}],result,"RawJSON"];
 Print["PROFILE PASS ",kind," ",order," ",ToString[value,InputForm]];Exit[0],fail["aborted"]];
