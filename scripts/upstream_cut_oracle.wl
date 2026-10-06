(* Cut-specific driver; upstream_cut_oracle.py reuses the bounded original
   runner, exact computational-source hashes and one-thread Wolfram launcher. *)
$HistoryLength = 0;
fail[message_] := (Print["ORACLE FAILURE: ", message]; Exit[1]);
assert[condition_, message_] := If[!TrueQ[condition], fail[message]];
base = StringReplace[DirectoryName[$InputFileName], RegularExpression["/+$"] -> ""];
root = Environment["AMFLOW_ROOT"];
run = Environment["AMFLOW_ORACLE_WORKDIR"];
oracleCase = Environment["AMFLOW_ORACLE_CASE"];
mode = Environment["AMFLOW_ORACLE_MODE"];
assert[StringQ[run] && DirectoryQ[run] && FileNames["*", run] === {}, "fresh oracle directory required"];
assert[StringStartsQ[run, base <> "/"], "work directory must be inside the oracle directory"];
assert[MemberQ[{"sample", "laurent"}, mode], "unknown oracle mode"];
assert[IntegerString[FileHash[FileNameJoin[{root, "AMFlow.m"}], "SHA256"], 16, 64] ===
 "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184", "upstream source mismatch"];
SetDirectory[run];
CheckAbort[
 Get[FileNameJoin[{root, "AMFlow.m"}]];
 AMFlow`Private`$WolframPath = Environment["AMFLOW_CONTROLLED_KERNEL"];
 AMFlow`SetReductionOptions["IBPReducer" -> "Kira"];
 AMFlow`SetAMFOptions["UseCache" -> False, "DESolver" -> "MMA", "RecursionMode" -> "AMF",
  "AMFMode" -> {"All"}, "D0" -> 4,
  "CacheName" -> FileNameJoin[{StringDrop[run, StringLength[base]+1], "cache"}]];
 CloseKernels[];
 AMFlow`AMFlowInfo["Family"] = cutoracle;
 AMFlow`AMFlowInfo["NThread"] = 1;
 AMFlow`AMFlowInfo["Conservation"] = {};
 AMFlow`AMFlowInfo["Loop"] = {r, k};
 AMFlow`AMFlowInfo["Leg"] = {p};
 sample = 1/13; dimension = 4-2eps;
 Switch[oracleCase,
  "mixed_cut_bubble",
   AMFlow`AMFlowInfo["Replacement"] = {p^2 -> 25};
   AMFlow`AMFlowInfo["Propagator"] = {r^2-1, (p-r)^2-4, k^2, (k-r)^2, (k-p)^2};
   AMFlow`AMFlowInfo["Prescription"] = {0, 1};
   AMFlow`AMFlowInfo["Cut"] = {1, 1, 0, 0, 0};
   targets = {j[cutoracle,1,1,1,1,0], j[cutoracle,2,1,1,1,0], j[cutoracle,1,2,1,1,0],
    j[cutoracle,1,1,1,1,-1], j[cutoracle,1,1,1,1,-2]};
   volume = Sqrt[6]/(25Pi) (25Pi/24)^eps Gamma[3/2]/Gamma[3/2-eps];
   bubble = Gamma[eps] Gamma[1-eps]^2/Gamma[2-2eps] Exp[I Pi eps];
   analytic = volume bubble {1, -7(1-2eps)/96-eps, -11(1-2eps)/192, 14,
     (121dimension-25)/(dimension-1)+75},
  "weighted_three_body",
   AMFlow`AMFlowInfo["Replacement"] = {p^2 -> 1};
   AMFlow`AMFlowInfo["Propagator"] = {r^2, k^2, (p-r-k)^2, (r+k)^2-2, (p-r)^2};
   AMFlow`AMFlowInfo["Prescription"] = {0, 0};
   AMFlow`AMFlowInfo["Cut"] = {1, 1, 1, 0, 0};
   targets = {j[cutoracle,1,1,1,1,0], j[cutoracle,1,1,1,2,0]};
   a = 1-eps;
   phi2 = (4Pi)^eps Gamma[a]/(8Pi Gamma[2a]);
   volume = phi2^2/(2Pi) Gamma[a] Gamma[2a]/Gamma[3a];
   analytic = volume {-Hypergeometric2F1[1,a,3a,1/2]/2, Hypergeometric2F1[2,a,3a,1/2]/4},
  _, fail["unknown cut case"]
 ];
 AMFlow`AMFlowInfo["Numeric"] = If[mode === "sample", {eps -> sample}, {}];
 requestedDigits = 24; expansionOrder = If[mode === "sample", 0, 4];
 config = AMFlow`GenerateNumericalConfig[requestedDigits, expansionOrder];
 measured = AbsoluteTiming[AMFlow`SolveIntegrals[targets, requestedDigits, expansionOrder]];
 assert[MatchQ[measured, {_?NumericQ, {__Rule}}], "malformed original result"];
 Put[measured[[2]], "oracle-result.wl"];
 assert[ContainsAll[Keys[measured[[2]]], targets], "missing cut targets"];
 answers = targets /. measured[[2]];
 powers = If[mode === "sample", {0}, Range[-4,0]];
 coefficients = If[mode === "sample", List /@ answers,
   Table[Coefficient[Expand[answer eps^4], eps, power+4], {answer, answers}, {power, powers}]];
 expected = If[mode === "sample", List /@ N[analytic /. eps -> sample,80],
   Table[N[SeriesCoefficient[expression,{eps,0,power}],80], {expression,analytic}, {power,powers}]];
 errors = MapThread[Abs[SetPrecision[#1,80]-#2]/Max[1,Abs[#2]]&, {coefficients,expected}, 2];
 assert[AllTrue[Flatten[errors],NumberQ], "nonnumeric cut comparison"];
 passed = TrueQ[Max[errors] < 10^-20];
 Export["oracle-metadata.json", <|
  "upstream_commit" -> "26005517a288086c4cb4d1b26d829691bc088485", "case" -> oracleCase,
  "mode" -> mode, "wolfram_version" -> $Version, "requested_digits" -> requestedDigits,
  "asserted_digits" -> 20, "targets" -> (ToString[#,InputForm]& /@ targets),
  "epsilon_samples_exact" -> (ToString[#,InputForm]& /@ config[[1]]),
  "working_decimal_precision" -> config[[2]], "x_order" -> config[[3]], "extra_x_order" -> config[[4]],
  "coefficient_powers" -> powers, "coefficients" -> Map[ToString[#,InputForm]&,coefficients,{2}],
  "expected" -> Map[ToString[#,InputForm]&,expected,{2}],
  "scaled_absolute_errors" -> Map[ToString[#,InputForm]&,errors,{2}],
  "solve_wall_seconds" -> measured[[1]], "analytic_check_passed" -> passed,
  "normalization" -> "Lorentz-invariant positive-energy phase space times d^D k/(i Pi^(D/2)) for each virtual loop; no EulerGamma factor, flux or particle symmetry factor",
  "normalization_conversion_to_rust" -> "1. Original guarded phase-volume ending uses 2*(-1)^(Lcut+1)*(4 Pi)^(-Lcut D/2) Im(I_uncut). This is used only for complete all-cut boundary components, never for the mixed graph.",
  "orientation" -> "Cut momenta r and p-r in mixed case; r,k,p-r-k in three-body case; p declared future timelike",
  "numerator_conversion" -> "Original extra slot is (k-p)^2 in mixed case, whereas native tests use k.p; scaleless virtual pinches imply first and second moments 14 I and [(121D-25)/(D-1)+75] I",
  "comparison" -> "Stored returned numbers padded only for residual computation; no additional precision claimed"
 |>, "RawJSON"];
 assert[passed, "cut oracle failed twenty-digit analytic agreement"];
 Print["ORACLE PASS: ", oracleCase, "; max scaled error = ", ToString[Max[errors],InputForm]];
 Exit[0], fail["original AMFlow aborted"]
];
fail["oracle ended without success"];
