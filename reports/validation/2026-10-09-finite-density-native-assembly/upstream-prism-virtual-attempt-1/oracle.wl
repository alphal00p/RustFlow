(* Independent validation-only two-cut prism virtual coefficients. No density
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
 sample=1/5;
 AMFlow`AMFlowInfo["Numeric"]={eps->sample};
 vars={d0,d1,d2,d3,d4,d5};
 n[s_]=((d0-d4-s)^2-(d0+d1-d5)*(1+s))/4;
 n0=n[0];
 np=(d0-d4)*(d4-d2-1)/2-(d0+d1-d5)/4;
 assert[Expand[(D[n[s],s]/.s->0)+D[n0,d4]*(d2-d4)-np]===0,"fixed-original-numerator derivative"];
 lower[polynomial_,powers_]:=Total[(Last[#]*Apply[j,Prepend[
    Join[Take[powers,6]-First[#],{0}],prismVirtual]])& /@ CoefficientRules[Expand[polynomial],vars]];
 expression0=lower[n0,{1,1,1,1,1,1,0}];
 expression1=lower[d4*np+n0*(d4-d2),{1,1,1,1,2,1,0}];
 targets=DeleteDuplicates[Cases[{expression0,expression1},j[__],Infinity]];
 assert[Length[targets]>0,"empty lowered target list"];
 Put[<|"numerator"->n0,"numerator_derivative"->np,"C0"->expression0,"C1"->expression1,"targets"->targets|>,"exact-targets.wl"];
 Export["definition.json",<|
   "role"->"independent validation-only virtual subgraph coefficients; not a complete density reference",
   "physical_input"->"triangular_prism; occupied pair [1,5]; K=P1,L=P3,q=P2,p=P4",
   "invariants_native_minkowski"->{"p^2=-1","q^2=0","p.q=-1/2"},
   "off_shell_parameter"->"s=q_E^2; (q-p)_E^2=0; derivative dq/ds=p-q at h=p_E^2=1",
   "virtual_jet"->"I_E(h,s)=h^(D-4)*(C0+C1*s/h+...); jet only, nonanalytic branches require common continuation",
   "raw_native_normalization"->"two loops d^D k/(i Pi^(D/2)); six physical denominators; no exp(EulerGamma eps)",
   "raw_euclidean_conversion"->"Cj_E=(4 Pi)^(-D) Cj_native",
   "expression_C0"->ToString[expression0,InputForm],"expression_C1"->ToString[expression1,InputForm],
   "target_count"->Length[targets],"epsilon"->ToString[sample,InputForm],"dimension"->ToString[4-2sample,InputForm],
   "common_thermal_endpoint_proved"->False,"complete_prism_reference"->False|>,"RawJSON"];
 profiles={}; values={};
 Do[
   AMFlow`SetAMFOptions["UseCache"->False,"DESolver"->"MMA","RecursionMode"->"AMF","D0"->4,
     "CacheName"->FileNameJoin[{"work","cache-"<>ToString[digits]}]];
   config=AMFlow`GenerateNumericalConfig[digits,0];
   timed=AbsoluteTiming[AMFlow`SolveIntegrals[targets,digits,0]];
   assert[MatchQ[timed,{_?NumericQ,{__Rule}}],"malformed automatic result"];
   assert[ContainsAll[Keys[timed[[2]]],targets],"missing requested target"];
   Put[timed[[2]],"result-"<>ToString[digits]<>".wl"];
   answer=({expression0,expression1}/.timed[[2]])/.eps->sample;
   assert[AllTrue[answer,NumberQ],"unresolved virtual coefficient"];
   AppendTo[values,answer];
   AppendTo[profiles,<|"requested_digits"->digits,"epsilon"->ToString[sample,InputForm],
     "native_values"->(ToString[#,InputForm]&/@answer),
     "euclidean_values"->(ToString[#,InputForm]&/@((4Pi)^(-(4-2sample))*answer)),
     "returned_precision"->(ToString[Precision[#],InputForm]&/@answer),
     "solve_wall_seconds"->timed[[1]],"working_precision"->config[[2]],
     "series_order"->config[[3]],"extra_series_order"->config[[4]]|>];
   Export["profiles.json",profiles,"RawJSON"];
   Print["PRISM VIRTUAL PROFILE COMPLETE: ",digits],
   {digits,{18,26}}];
 residuals=MapThread[Abs[SetPrecision[#1,80]-SetPrecision[#2,80]]/Max[1,Abs[SetPrecision[#2,80]]]&,
   {values[[1]],values[[2]]}];
 passed=TrueQ[Max[residuals]<10^-12];
 Export["comparison.json",<|"kind"->"independent upstream virtual coefficient precision refinement",
   "status"->If[passed,"passed","failed"],"scaled_absolute_threshold"->"1e-12",
   "residuals"->(ToString[#,InputForm]&/@residuals),"profiles"->profiles,
   "wolfram_version"->$Version,"upstream_commit"->"26005517a288086c4cb4d1b26d829691bc088485",
   "complete_prism_reference"->False,"numerical_oracle_records_compared"->0,
   "scope"->"C0 and formal off-shell jet C1 only; no full cut assembly, Laurent coefficients, or thermal/massless endpoint proof"|>,"RawJSON"];
 assert[passed,"precision refinement below 12 scaled absolute digits"];
 Print["PRISM VIRTUAL REFINEMENT PASS"];
 Exit[0],fail["upstream calculation aborted"]];
fail["driver ended without success"];
