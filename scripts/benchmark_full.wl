(* Optional full-workflow timing. Requires Wolfram plus a configured reducer.
   Run in a NEW working directory; AMFLOW_ROOT points at the pinned checkout.
   Example: AMFLOW_ROOT=/path/to/amflow AMFLOW_BENCH_CASE=paper \
            AMFLOW_IBP_REDUCER=Blade wolframscript -file benchmark_full.wl
   This harness has not been executed in the current environment. *)
root = Environment["AMFLOW_ROOT"];
If[!StringQ[root] || !FileExistsQ[FileNameJoin[{root, "AMFlow.m"}]],
  Print["AMFLOW_ROOT must name the pinned AMFlow checkout."]; Exit[2]];
If[IntegerString[FileHash[FileNameJoin[{root,"AMFlow.m"}],"SHA256"],16,64] =!=
   "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184",
  Print["AMFlow.m differs from the pinned upstream commit."]; Exit[2]];
case = Environment["AMFLOW_BENCH_CASE"];
If[!MemberQ[{"bubble", "sunset", "paper"}, case],
  Print["AMFLOW_BENCH_CASE must be bubble, sunset, or paper."]; Exit[2]];
reducer = Environment["AMFLOW_IBP_REDUCER"];
If[!MemberQ[{"Blade", "FiniteFlow+LiteRed", "FIRE+LiteRed", "Kira"}, reducer],
  Print["AMFLOW_IBP_REDUCER must name a configured upstream reducer."]; Exit[2]];
If[FileNames["amflow-benchmark-*"] =!= {},
  Print["Use a new working directory for each run."]; Exit[2]];
Get[FileNameJoin[{root, "AMFlow.m"}]];
SetReductionOptions["IBPReducer" -> reducer];
SetAMFOptions["UseCache" -> False, "RecursionMode" -> "AMF"];
AMFlowInfo["Family"] = bench;
AMFlowInfo["NThread"] = 1;
AMFlowInfo["Conservation"] = {};
Switch[case,
 "bubble",
 AMFlowInfo["Loop"] = {l}; AMFlowInfo["Leg"] = {p};
 AMFlowInfo["Replacement"] = {p^2 -> 1};
 AMFlowInfo["Propagator"] = {l^2, (l+p)^2};
 AMFlowInfo["Numeric"] = {}; target = {j[bench,1,1]}; order = 2,
 "sunset",
 AMFlowInfo["Loop"] = {l1,l2}; AMFlowInfo["Leg"] = {};
 AMFlowInfo["Replacement"] = {};
 AMFlowInfo["Propagator"] = {l1^2-1,l2^2-1,(l1+l2)^2};
 AMFlowInfo["Numeric"] = {}; target = {j[bench,1,1,1]}; order = 4,
 "paper",
 AMFlowInfo["Loop"] = {l1,l2}; AMFlowInfo["Leg"] = {p1,p2,p3,p4};
 AMFlowInfo["Conservation"] = {p4 -> -p1-p2-p3};
 AMFlowInfo["Replacement"] = {p1^2->0,p2^2->0,p3^2->msq,p4^2->msq,
   (p1+p2)^2->s,(p1+p3)^2->t};
 AMFlowInfo["Propagator"] = {l1^2,(l1+p1)^2,(l1+p1+p2)^2,l2^2,
   (l2+p3)^2-msq,(l2+p3+p4)^2,(l1+l2)^2,(l1-p3)^2,(l2+p1)^2};
 AMFlowInfo["Numeric"] = {s->30,t->-10/3,msq->1};
 target = Table[j[bench,1,1,1,1,1,1,1,-3+i,-i],{i,0,3}]; order = 4
];
cpuStart = TimeUsed[];
measured = CheckAbort[AbsoluteTiming[SolveIntegrals[target,20,order]], $Aborted];
If[measured === $Aborted || !MatchQ[measured, {_?NumericQ, {__Rule}}],
  Print["Automatic evaluation failed; no successful timing recorded."]; Exit[1]];
Put[measured[[2]], "amflow-benchmark-result.wl"];
Export["amflow-benchmark-metadata.json", <|
 "upstream_commit"->"26005517a288086c4cb4d1b26d829691bc088485",
 "case"->case,"wolfram_version"->$Version,"system"->$System,
 "reducer"->reducer,"threads"->1,"requested_digits"->20,"epsilon_through"->0,
 "solve_wall_seconds"->measured[[1]],"main_kernel_cpu_seconds"->TimeUsed[]-cpuStart,
 "cache"->False,"verified_accuracy"->False,
 "note"->"Compare coefficients and independently refine before interpreting timings."|>, "RawJSON"];
Exit[0];
