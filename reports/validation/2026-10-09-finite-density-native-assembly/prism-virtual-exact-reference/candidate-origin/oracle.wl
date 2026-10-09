(* Independent validation-only all two-cut prism virtual coefficients. No density
   amplitude, production boundary, oracle answer or thermal prescription. *)
$HistoryLength=0;
fail[m_]:=(Print["PRISM VIRTUAL FAILURE: ",m];Exit[1]);
assert[c_,m_]:=If[!TrueQ[c],fail[m]];
base=DirectoryName[$InputFileName];
root=Environment["AMFLOW_ROOT"];
work=Environment["AMFLOW_ORACLE_WORKDIR"];
assert[DirectoryQ[root]&&DirectoryQ[work],"isolated directories missing"];
assert[FileNames["*",work]==={},"working directory not fresh"];
assert[IntegerString[FileHash[FileNameJoin[{root,"AMFlow.m"}],"SHA256"],16,64]===
 "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184","upstream hash"];
SetDirectory[work];
CheckAbort[
 Get[FileNameJoin[{root,"AMFlow.m"}]];
 AMFlow`Private`$WolframPath=Environment["AMFLOW_CONTROLLED_KERNEL"];
 AMFlow`SetReductionOptions["IBPReducer"->"Kira"];
 CloseKernels[];
 AMFlow`AMFlowInfo["Family"]=prismVirtual;
 AMFlow`AMFlowInfo["NThread"]=1;
 AMFlow`AMFlowInfo["Loop"]={k,l};
 AMFlow`AMFlowInfo["Leg"]={p,q};
 AMFlow`AMFlowInfo["Conservation"]={};
 AMFlow`AMFlowInfo["Replacement"]={p^2->-1,q^2->0,p*q->-1/2};
 AMFlow`AMFlowInfo["Propagator"]={k^2,l^2,(k-p)^2,(l-p)^2,(k-q)^2,(k-l)^2,(l-q)^2};
 sample=eps;
 AMFlow`AMFlowInfo["Numeric"]={};
 vars={d0,d1,d2,d3,d4,d5};
 n[s_]=((d0-d4-s)^2-(d0+d1-d5)*(1+s))/4;
 n0=n[0];
 np=(d0-d4)*(d4-d2-1)/2-(d0+d1-d5)/4;
 assert[Expand[(D[n[s],s]/.s->0)+D[n0,d4]*(d2-d4)-np]===0,"fixed-original-numerator derivative"];
 lower[polynomial_,powers_]:=Total[(Last[#]*Apply[j,Prepend[
    Join[Take[powers,6]-First[#],{0}],prismVirtual]])& /@ CoefficientRules[Expand[polynomial],vars]];
 expression0=lower[n0,{1,1,1,1,1,1,0}];
 expression1=lower[d4*np+n0*(d4-d2),{1,1,1,1,2,1,0}];
 nb0=(1+(d1-1-d3)*(d0-d4))/4;
 nbp=1/2+(d1-1-d3)*(d4-d2-1)/4;
 nc=((d0+d4)^2+(d0-d1+d5)*(d4+d2))/4;
 expressionB0=lower[nb0,{1,1,1,1,1,1,0}];
 expressionB1=lower[d4*nbp+nb0*(d4-d2),{1,1,1,1,2,1,0}];
 expressionC=lower[nc,{1,1,1,1,2,1,0}];
 expressions={expression0,expression1,expressionB0,expressionB1,expressionC};
 labels={"C0A","C1A","C0B","C1B","CC"};
 targets=DeleteDuplicates[Cases[expressions,j[__],Infinity]];
 assert[Length[targets]>0,"empty lowered target list"];
 Put[<|"numerator"->n0,"numerator_derivative"->np,"C0"->expression0,"C1"->expression1,"targets"->targets|>,"exact-targets.wl"];
 Export["definition.json",<|
   "role"->"independent validation-only virtual subgraph coefficients; not a complete density reference",
   "physical_input"->"triangular_prism; pair A [1,5], pair B [1,7], pair C [5,7]; exact maps in finite-density-prism-reference-construction.md",
   "invariants_native_minkowski"->{"p^2=-1","q^2=0","p.q=-1/2"},
   "off_shell_parameter"->"s=q_E^2; (q-p)_E^2=0; derivative dq/ds=p-q at h=p_E^2=1",
   "virtual_jet"->"I_E(h,s)=h^(D-4)*(C0+C1*s/h+...); jet only, nonanalytic branches require common continuation",
   "raw_native_normalization"->"two loops d^D k/(i Pi^(D/2)); six physical denominators; no exp(EulerGamma eps)",
   "raw_euclidean_conversion"->"(4 Pi)^(-D) times native [C0A,C1A,C0B,C1B,-CC]",
   "coefficient_labels"->labels,"exact_expressions"->(ToString[#,InputForm]&/@expressions),
   "target_count"->Length[targets],"epsilon"->ToString[sample,InputForm],"dimension"->ToString[4-2sample,InputForm],
   "common_thermal_endpoint_proved"->False,"complete_prism_reference"->False|>,"RawJSON"];
 AMFlow`SetAMFOptions["UseCache"->False,"CacheName"->FileNameJoin[{"work","combination-reduction"}]];
 reductionTiming=AbsoluteTiming[AMFlow`BlackBoxReduce[targets,{}]];
 red=reductionTiming[[2]];
 assert[MatchQ[red,{_List,_List}],"malformed exact numerator reduction"];
 Put[red,"exact-combination-reduction.wl"];
 rules=(First[#]->Together[Last[#].red[[1]]])&/@red[[2]];
 reducedExpressions=Factor[Together[#]]&/@(expressions/.rules);
 targets=DeleteDuplicates[Cases[reducedExpressions,j[__],Infinity]];
 Put[<|"labels"->labels,"original"->expressions,"reduced"->reducedExpressions,"surviving_masters"->targets|>,"reduced-combinations.wl"];
 Print["PRISM COMPLETE COMBINATIONS REDUCED: ",Length[red[[1]]]," masters before coefficient cancellation; ",Length[targets]," surviving."];

 expectedMasters={j[prismVirtual,0,1,0,1,1,1,0],j[prismVirtual,0,1,1,0,0,1,0],j[prismVirtual,1,1,1,1,0,0,0]};
 assert[Sort[targets]===Sort[expectedMasters],"generic regulator requires additional virtual master; preserve reduction and stop"];
 bubble=Gamma[eps] Gamma[1-eps]^2/Gamma[2-2eps];
 masterGamma={bubble Gamma[2eps] Gamma[1-2eps]^2/Gamma[2-3eps],-Gamma[1-eps]^3 Gamma[2eps-1]/Gamma[3-3eps],bubble^2};
 coefficients=reducedExpressions/.Thread[expectedMasters->masterGamma];
 Put[<|"labels"->labels,"masters"->expectedMasters,"master_gamma"->masterGamma,"native_coefficients"->coefficients|>,"symbolic-gamma-coefficients.wl"];
 Export["symbolic-reduction-status.json",<|"status"->"symbolic_virtual_coefficients_derived","epsilon"->"symbolic","reduction_wall_seconds"->reductionTiming[[1]],"surviving_master_count"->Length[targets],"complete_prism_reference"->False,"thermal_mass_derivative_continuation_proved"->False,"native_predictions_read"->0,"oracle_records_read"->0,"coefficient_labels"->labels,"native_coefficients"->(ToString[#,InputForm]&/@coefficients)|>,"RawJSON"];
 Print["PRISM SYMBOLIC VIRTUAL REDUCTION PASS"];
 Exit[0],fail["upstream calculation aborted"]];
fail["driver ended without success"];
