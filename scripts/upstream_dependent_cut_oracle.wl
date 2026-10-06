(* Uses unchanged pinned AMFlow computational sources and the bounded runner. *)
$HistoryLength = 0;
fail[message_] := (Print["ORACLE FAILURE: ", message]; Exit[1]);
assert[condition_, message_] := If[!TrueQ[condition], fail[message]];
base = StringReplace[DirectoryName[$InputFileName], RegularExpression["/+$"] -> ""];
root = Environment["AMFLOW_ROOT"];
run = Environment["AMFLOW_ORACLE_WORKDIR"];
oracleCase = Environment["AMFLOW_ORACLE_CASE"];
assert[Environment["AMFLOW_ORACLE_MODE"] === "sample", "this oracle compares finite-epsilon samples"];
assert[DirectoryQ[run] && FileNames["*", run] === {}, "fresh directory required"];
assert[StringStartsQ[run, base <> "/"], "work directory must belong to the runner"];
assert[IntegerString[FileHash[FileNameJoin[{root,"AMFlow.m"}],"SHA256"],16,64] ===
 "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184", "source mismatch"];
SetDirectory[run];
CheckAbort[
 Get[FileNameJoin[{root,"AMFlow.m"}]];
 AMFlow`Private`$WolframPath = Environment["AMFLOW_CONTROLLED_KERNEL"];
 AMFlow`SetReductionOptions["IBPReducer" -> "Kira"];
 firstMass = Switch[oracleCase,"dependent_child_two",2,"dependent_child_three",3,_,fail["unknown case"]];
 placement = "All";
 AMFlow`SetAMFOptions["UseCache" -> False,"DESolver" -> "MMA","RecursionMode" -> "AMF",
  "AMFMode" -> {placement},"D0" -> 4,
  "CacheName" -> FileNameJoin[{StringDrop[run,StringLength[base]+1],"cache"}]];
 CloseKernels[];
 AMFlow`AMFlowInfo["Family"] = partialcut;
 AMFlow`AMFlowInfo["NThread"] = 1;
 AMFlow`AMFlowInfo["Conservation"] = {};
 AMFlow`AMFlowInfo["Loop"] = {r,l,k};
 AMFlow`AMFlowInfo["Leg"] = {p};
 AMFlow`AMFlowInfo["Replacement"] = {p^2 -> 1};
 AMFlow`AMFlowInfo["Propagator"] = {r^2,l^2,(p-r-l)^2,k^2-firstMass,(k-r-l)^2-5,
  (p-r)^2,(p-l)^2,(k-r)^2,(p-k)^2};
 AMFlow`AMFlowInfo["Prescription"] = {0,0,1};
 AMFlow`AMFlowInfo["Cut"] = {1,1,1,0,0,0,0,0,0};
 sample = 1/13;
 AMFlow`AMFlowInfo["Numeric"] = {eps -> sample};
 targets = {j[partialcut,1,1,1,1,1,0,0,0,0]};
 requestedDigits = 24;
 measured = AbsoluteTiming[AMFlow`SolveIntegrals[targets,requestedDigits,0]];
 assert[MatchQ[measured,{_?NumericQ,{__Rule}}],"malformed result"];
 Put[measured[[2]],"oracle-result.wl"];
 answers = targets /. measured[[2]];
 a = 1-sample;
 phi2 = (4Pi)^sample Gamma[a]/(8Pi Gamma[2a]);
 volume = phi2^2/(2Pi) Gamma[a] Gamma[2a]/Gamma[3a];
 terms = Table[Pochhammer[sample,n]/(n! 5^n) Sum[
  Binomial[n,j] (5-firstMass)^(n-j) n! j!/(n+j+1)! Pochhammer[a,j]/Pochhammer[3a,j],{j,0,n}],{n,0,219}];
 expected = {N[volume Gamma[sample] 5^-sample Total[terms],80]};
 tailBound = N[Abs[volume Gamma[sample] 5^-sample] (3/5)^220/(1-3/5),80];
 errors = MapThread[Abs[SetPrecision[#1,80]-#2]/Max[1,Abs[#2]]&,{answers,expected}];
 passed = TrueQ[Max[errors]+tailBound < 10^-20];
 configs = Table[
  cfg = Get[file];
  props = "Propagator" /. cfg;
  assert[ListQ[props],"invalid recorded AMFlow propagators"];
  <|"path" -> StringDrop[file,StringLength[run]+1],
    "loops" -> ToString["Loop" /. cfg,InputForm],
    "propagators" -> ToString[props,InputForm],
    "eta_positions_one_based" -> Flatten[Position[(!FreeQ[#,AMFlow`Private`$Eta]& /@ props),True]]|>,
  {file,Sort[Select[FileNames["config",run,Infinity],FileType[#] === File&]]}];
 assert[Length[configs]>0,"missing actual upstream placement metadata"];
 assert[AnyTrue[configs,#["eta_positions_one_based"] === If[placement==="All",{4,5},{4}]&],"top-system placement differs from the intended comparison"];
 Export["oracle-metadata.json",<|
  "upstream_commit" -> "26005517a288086c4cb4d1b26d829691bc088485",
  "case" -> oracleCase,"mode" -> "sample","requested_amf_mode" -> placement,
  "actual_system_configurations" -> configs,"epsilon_exact" -> ToString[sample,InputForm],
  "requested_digits" -> requestedDigits,"asserted_digits" -> 20,"wolfram_version" -> $Version,
  "coefficients" -> (ToString[#,InputForm]& /@ answers),
  "expected" -> (ToString[#,InputForm]& /@ expected),
  "scaled_absolute_errors" -> (ToString[#,InputForm]& /@ errors),
  "independent_series_absolute_tail_bound" -> ToString[tailBound,InputForm],
  "solve_wall_seconds" -> measured[[1]],"analytic_check_passed" -> passed,
  "normalization" -> "Complete massless three-body Lorentz-invariant phase space times d^D k/(i Pi^(D/2)); no EulerGamma, flux or symmetry factor",
  "independence" -> "Comparison-only double-beta/Feynman-parameter series; no numerical reference boundary enters the evaluator",
  "comparison_scope" -> "Independent child of native overcomplete inventory. Combined value is child(mass=3)-child(mass=2). No claim that upstream accepts the original dependent inventory.",
  "first_squared_mass" -> firstMass
 |>,"RawJSON"];
 assert[passed,"twenty-digit nonfactorized comparison failed"];
 Print["ORACLE PASS: ",oracleCase,"; scaled error ",Max[errors]];
 Exit[0],fail["original AMFlow aborted"]
];
fail["oracle ended without success"];
