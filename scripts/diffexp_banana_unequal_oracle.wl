(* Independently written driver for the unchanged, pinned DiffExp banana.
   DIFFEXP_ROOT: pinned checkout; DIFFEXP_OUTPUT: NEW output directory.
   DIFFEXP_BANANA_ROUTE: shortcut (default) or full. Both initialize from
   analytic gamma factors, never from numerical reference coefficients.
   Use scripts/run_diffexp_banana_unequal.py to enforce time/RSS bounds.
*)
$HistoryLength=0;
driverFile=ExpandFileName[$InputFileName];
fail[m_]:=(Print["ORACLE FAILURE: ",m];Exit[1]);
requiredEnvironment[name_]:=Module[{value=Environment[name]},
 If[!StringQ[value]||StringLength[value]==0,fail["Set "<>name<>" before running"]];
 ExpandFileName[value]];
base=requiredEnvironment["DIFFEXP_OUTPUT"];root=requiredEnvironment["DIFFEXP_ROOT"];
route=Environment["DIFFEXP_BANANA_ROUTE"];
If[!StringQ[route]||StringLength[route]==0,route="shortcut"];
If[!MemberQ[{"full","shortcut"},route],fail["DIFFEXP_BANANA_ROUTE must be full or shortcut"]];
shortcut=route==="shortcut";
If[FileExistsQ[base]||DirectoryQ[base],fail["DIFFEXP_OUTPUT must be new"]];
If[CreateDirectory[base,CreateIntermediateDirectories->True]===$Failed,fail["Cannot create output"]];
repository=DirectoryName[DirectoryName[driverFile]];
fixture=Import[FileNameJoin[{repository,"fixtures","diffexp","banana-unequal.json"}],"RawJSON"];
equalFixture=Import[FileNameJoin[{repository,"fixtures","diffexp","banana-equal.json"}],"RawJSON"];
If[!AssociationQ[fixture]||!AssociationQ[equalFixture],fail["Missing banana fixtures"]];

sha[f_]:=IntegerString[FileHash[f,"SHA256"],16,64];str[v_]:=ToString[v,InputForm];
scalar[v_]:=<|"real_wolfram"->str[Re[v]],"imaginary_wolfram"->str[Im[v]],"precision_wolfram"->str[Precision[v]],"accuracy_wolfram"->str[Accuracy[v]]|>;
writeJSON[path_,data_]:=If[Export[path,data,"RawJSON"]===$Failed,fail["JSON export failed"]];
If[sha[FileNameJoin[{root,"DiffExp.m"}]]=!="67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7",fail["Pinned package mismatch"]];
Do[If[!FileExistsQ[FileNameJoin[{root,entry["file"]}]]||sha[FileNameJoin[{root,entry["file"]}]]=!=entry["sha256"],fail["Matrix hash mismatch"]],{entry,Join[fixture["source_files"],equalFixture["matrices"]]}];
Get[FileNameJoin[{root,"DiffExp.m"}]];
rawBoundary={"?","?",\[Epsilon]*(1+3*\[Epsilon])*(1+4*\[Epsilon])*E^(3*EulerGamma*\[Epsilon])*(-4*Gamma[\[Epsilon]]^3/t+6*(-1/t)^(1+\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^2*Gamma[\[Epsilon]]^3/Gamma[-2*\[Epsilon]]+8*(-1/t)^(1+2*\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^3*Gamma[\[Epsilon]]*Gamma[2*\[Epsilon]]/Gamma[-3*\[Epsilon]]+3*(-1/t)^(1+3*\[Epsilon])*\[Epsilon]*Gamma[-\[Epsilon]]^4*Gamma[3*\[Epsilon]]/Gamma[-4*\[Epsilon]]),E^(3*EulerGamma*\[Epsilon])*\[Epsilon]^3*Gamma[\[Epsilon]]^3};
records=<||>;
record[name_,timed_,expected_]:=Module[{value=timed[[2]]},
 If[Dimensions[value[[2]]]=!=expected||!AllTrue[Flatten[value[[2]]],NumberQ],fail["Malformed stage "<>name]];
 Put[value,FileNameJoin[{base,name<>".wl"}]];
 AssociateTo[records,name-><|"seconds"->timed[[1]],"point"->str[value[[1]]],"values"->Map[scalar,value[[2]],{2}],"reported_errors"->If[Length[value]>=3,Map[str,value[[3]],{2}],Null]|>];
 writeJSON[FileNameJoin[{base,"stages.json"}],records];Print["STAGE DONE: ",name," in ",timed[[1]]," seconds"]];
CheckAbort[
 DiffExp`LoadConfiguration[{DiffExp`MatrixDirectory->FileNameJoin[{root,"Examples","Banana_EqualMass_Matrices"}],DiffExp`EpsilonOrder->4,WorkingPrecision->500,DiffExp`ChopPrecision->250,DiffExp`ExpansionOrder->If[shortcut,80,70],DiffExp`DivisionOrder->If[shortcut,2,4],DiffExp`RadiusOfConvergence->1,DiffExp`UseMobius->True,DiffExp`UsePade->True,DiffExp`DeltaPrescriptions->{t-16+I*\[Delta]},DiffExp`Verbosity->1,"Parallel"->False}];
 prepared=AbsoluteTiming[DiffExp`PrepareBoundaryConditions[rawBoundary,{t->-1/x}]];
 record["equal_minus1",AbsoluteTiming[DiffExp`TransportTo[prepared[[2]],{t->-1}]],{4,5}];
 equalMinus1=Get[FileNameJoin[{base,"equal_minus1.wl"}]];
 equalPoint=If[shortcut,50,1/2];equalName=If[shortcut,"equal_50","equal_half"];
 record[equalName,AbsoluteTiming[DiffExp`TransportTo[equalMinus1,{t->x},equalPoint]],{4,5}];
 equalEndpoint=Get[FileNameJoin[{base,equalName<>".wl"}]];
 mapping={1,1,1,1,1,1,2,2,2,2,3,4,4,4,4};
 bc={<|psq->equalPoint,mm1->1,mm2->1,mm3->1,mm4->1|>,equalEndpoint[[2,mapping]]};
 DiffExp`LoadConfiguration[{DiffExp`MatrixDirectory->FileNameJoin[{root,"Examples","Banana_Matrices"}],DiffExp`EpsilonOrder->4,WorkingPrecision->1000,DiffExp`ChopPrecision->500,DiffExp`ExpansionOrder->70,DiffExp`DivisionOrder->4,DiffExp`RadiusOfConvergence->10,DiffExp`UseMobius->True,DiffExp`UsePade->True,DiffExp`Verbosity->1,"Parallel"->False}];
 massName=If[shortcut,"mass_deformation_50","mass_deformation_half"];
 record[massName,AbsoluteTiming[DiffExp`TransportTo[bc,<|psq->equalPoint,mm1->1+x,mm2->1+x/2,mm3->1+x/3,mm4->1|>,1]],{15,5}];
 final=Get[FileNameJoin[{base,massName<>".wl"}]];
 savedMetadata=<||>;
 If[!shortcut,
  DiffExp`UpdateConfiguration[DiffExp`DeltaPrescriptions->{(Sqrt[2]+Sqrt[3/2]+Sqrt[4/3]+1)^2-psq-I*\[Delta]}];
  physical=AbsoluteTiming[DiffExp`TransportTo[final,{psq->x,mm1->2,mm2->3/2,mm3->4/3,mm4->1},50,True]];
  record["physical_50",{physical[[1]],physical[[2,1]]},{15,5}];
  Put[physical[[2]],FileNameJoin[{base,"physical_saved_series.wl"}]];
  saved=AbsoluteTiming[DiffExp`ToPiecewise[physical[[2]],True]];
  at10=AbsoluteTiming[Table[saved[[2,i,j]][10],{i,15},{j,5}]];
  If[!AllTrue[Flatten[at10[[2]]],NumberQ],fail["Saved series at10"]];
  savedMetadata=<|"saved_pade_construction_seconds"->saved[[1]],"saved_at10_evaluation_seconds"->at10[[1]],"saved_at10_values"->Map[scalar,at10[[2]],{2}]|>;
  final=Get[FileNameJoin[{base,"physical_50.wl"}]]];
 refs=fixture["reference"]["paper_B11"];
 exactDecimal[s_]:=Module[{parts=StringSplit[s,"."]},If[Length[parts]==1,ToExpression[s],ToExpression[StringJoin[parts]]/10^StringLength[parts[[2]]]]];
 expected=(exactDecimal[#["real"]]+I*exactDecimal[#["imaginary"]])&/@refs;
 paperDifference=Max[Abs[final[[2,11]]-expected]];
 (* Decimal paper values are comparison-only, never inputs to a boundary. *)
 If[!TrueQ[paperDifference<10^-20],fail[{"B11 paper comparison",paperDifference}]];
 writeJSON[FileNameJoin[{base,"result.json"}],<|"status"->"passed","case"->"unequal_mass_banana","route"->route,"upstream_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd","driver_sha256"->sha[driverFile],"package_sha256"->sha[FileNameJoin[{root,"DiffExp.m"}]],"wolfram_version"->$Version,"settings"-><|"working_precision"->1000,"chop_precision"->500,"order"->70,"division_order"->4,"radius_of_convergence"->10,"mobius"->True,"pade"->True,"epsilon_order"->4|>,"equal_seed_settings"-><|"working_precision"->500,"chop_precision"->250,"expansion_order"->If[shortcut,80,70],"division_order"->If[shortcut,2,4],"radius"->1|>,"preparation_seconds"->prepared[[1]],"stages"->records,"saved_series"->savedMetadata,"max_absolute_B11_paper_difference"->str[paperDifference],"B11_asserted_digits"->20,"scope"->"Complete15-master epsilon0..4 original transport, analytic boundary; endpoint B11 compared to paper data. All-master independent comparison is performed separately.","notebook_deviations"->{"Verbosity1 rather than2; explicit single-process settings; resource limits supplied by external runner. Full route uses notebook settings; shortcut is the separately recorded equal-t50 route."},"input_error_handling"->"Notebook's equal-to-unequal boundary duplicates only values and omits its numerical error array. Equal-stage reported errors are retained separately in this report."|>];
 Print["ORACLE PASS ",str[paperDifference]];Exit[0],fail["DiffExp aborted"]];
