(* Optional live oracle against pinned AMFlow 2.0. Run via upstream_oracle.py,
   which copies this script next to its fresh working directory. One AMFlow
   worker is used; the parent and solver kernels may coexist. *)
$HistoryLength = 0;
fail[message_, code_: 1] := (Print["ORACLE FAILURE: ", message]; Exit[code]);
assert[condition_, message_] := If[!TrueQ[condition], fail[message]];
base = StringReplace[DirectoryName[$InputFileName], RegularExpression["/+$"] -> ""];
root = Environment["AMFLOW_ROOT"];
assert[StringQ[root] && DirectoryQ[root], "AMFLOW_ROOT must name the configured upstream copy."];
run = Environment["AMFLOW_ORACLE_WORKDIR"];
oracleCase = Environment["AMFLOW_ORACLE_CASE"];
mode = Environment["AMFLOW_ORACLE_MODE"];
assert[StringQ[run] && DirectoryQ[run], "A fresh AMFLOW_ORACLE_WORKDIR is required."];
assert[StringStartsQ[run, base <> "/"], "Working directory must be inside the oracle directory."];
assert[FileNames["*", run] === {}, "Refusing a nonempty working directory."];
assert[MemberQ[{"bubble", "sunset"}, oracleCase], "Unknown case."];
assert[MemberQ[{"sample", "laurent"}, mode], "Unknown mode."];
assert[IntegerString[FileHash[FileNameJoin[{root, "AMFlow.m"}], "SHA256"], 16, 64] ===
  "76feadd3990586dbba22c96f515328fd79b64d6767ed021f2f80c9e4a2b98184",
  "AMFlow.m differs from the pinned original."];
SetDirectory[run];
CheckAbort[
  (* Initial FiniteFlow dependency checks only print status. The Kira adapter
     selected immediately below needs neither FiniteFlow nor LiteRed. *)
  Get[FileNameJoin[{root, "AMFlow.m"}]];
  (* Runtime launcher configuration; generated solver scripts stay unchanged. *)
  controlledKernel = Environment["AMFLOW_CONTROLLED_KERNEL"];
  assert[StringQ[controlledKernel] && FileExistsQ[controlledKernel],
    "Missing one-thread child kernel launcher."];
  AMFlow`Private`$WolframPath = controlledKernel;
  AMFlow`SetReductionOptions["IBPReducer" -> "Kira"];
  AMFlow`SetAMFOptions["UseCache" -> False, "DESolver" -> "MMA",
    "RecursionMode" -> "AMF", "D0" -> 4,
    (* Original AMFlow resolves CacheName relative to this script, not Directory[]. *)
    "CacheName" -> FileNameJoin[{StringDrop[run, StringLength[base]+1], "cache"}]];
  CloseKernels[];
  AMFlow`AMFlowInfo["Family"] = oracle;
  AMFlow`AMFlowInfo["NThread"] = 1;
  AMFlow`AMFlowInfo["Conservation"] = {};
  Switch[oracleCase,
    "bubble",
      AMFlow`AMFlowInfo["Loop"] = {l};
      AMFlow`AMFlowInfo["Leg"] = {p};
      AMFlow`AMFlowInfo["Replacement"] = {p^2 -> -1};
      AMFlow`AMFlowInfo["Propagator"] = {l^2, (l+p)^2};
      target = {j[oracle, 1, 1]}; loops = 1; sample = 1/10;
      analytic = Gamma[eps] Gamma[1-eps]^2/Gamma[2-2eps];
      expectedCoefficients = {0, 1, 2-EulerGamma},
    "sunset",
      AMFlow`AMFlowInfo["Loop"] = {l1, l2};
      AMFlow`AMFlowInfo["Leg"] = {};
      AMFlow`AMFlowInfo["Replacement"] = {};
      AMFlow`AMFlowInfo["Propagator"] = {l1^2-1, l2^2-1, (l1+l2)^2};
      target = {j[oracle, 1, 1, 1]}; loops = 2; sample = 3/4;
      analytic = Gamma[eps]^2/((1-eps) (1-2eps));
      expectedCoefficients = {0, 0, 1, 3-2EulerGamma,
        7-6EulerGamma+2EulerGamma^2+Pi^2/6}
  ];
  AMFlow`AMFlowInfo["Numeric"] = If[mode === "sample", {eps -> sample}, {}];
  requestedDigits = 24;
  (* Original API order is relative to the leading epsilon^(-2 loops). *)
  expansionOrder = If[mode === "sample", 0, 2loops];
  config = AMFlow`GenerateNumericalConfig[requestedDigits, expansionOrder];
  cpuStart = TimeUsed[];
  measured = AbsoluteTiming[
    AMFlow`SolveIntegrals[target, requestedDigits, expansionOrder]];
  assert[MatchQ[measured, {_?NumericQ, {__Rule}}], "Malformed automatic result."];
  Put[measured[[2]], "oracle-result.wl"];
  (* The original Kira adapter may also return solved master integrals. *)
  assert[ContainsAll[Keys[measured[[2]]], target], "A requested target is missing."];
  answer = target[[1]] /. measured[[2]];
  If[mode === "sample",
    coefficients = {answer};
    expected = {N[analytic /. eps -> sample, 80]};
    powers = {0},
    powers = Range[-2loops, 0];
    coefficients = Table[Coefficient[Expand[answer eps^(2loops)], eps, power+2loops],
      {power, powers}];
    expected = N[expectedCoefficients, 80]
  ];
  assert[AllTrue[coefficients, NumberQ], "Unresolved or nonnumeric coefficients."];
  (* Expose the stored number's residual without significance-arithmetic
     cancellation to an inexact zero. Padding adds no accuracy guarantee. *)
  errors = MapThread[Abs[SetPrecision[#1, 80]-#2]/Max[1, Abs[#2]]&,
    {coefficients, expected}];
  assert[AllTrue[errors, NumberQ], "Nonnumeric error estimates."];
  Put[<|"expected" -> expected, "scaled_absolute_errors" -> errors|>,
    "oracle-comparison.wl"];
  passed = TrueQ[Max[errors] < 10^-20];
  Export["oracle-metadata.json", <|
    "upstream_commit" -> "26005517a288086c4cb4d1b26d829691bc088485",
    "case" -> oracleCase, "mode" -> mode, "wolfram_version" -> $Version,
    "reducer" -> "Kira 3.1 dad964cd79a39a9fca728bb2606f59537a95109c",
    "desolver" -> "MMA", "workers" -> 1, "cache" -> False,
    "cache_directory" -> FileNameJoin[{run, "cache"}],
    "requested_digits" -> requestedDigits, "asserted_digits" -> 20,
    "epsilon_samples_exact" -> (ToString[#, InputForm]& /@ config[[1]]),
    "working_decimal_precision" -> config[[2]], "x_order" -> config[[3]],
    "extra_x_order" -> config[[4]], "leading_epsilon_power" -> config[[5]],
    "coefficient_powers" -> powers,
    "coefficients" -> (ToString[#, InputForm]& /@ coefficients),
    "expected" -> (ToString[#, InputForm]& /@ expected),
    "scaled_absolute_errors" -> (ToString[#, InputForm]& /@ errors),
    "solve_wall_seconds" -> measured[[1]],
    "parent_kernel_cpu_seconds" -> TimeUsed[]-cpuStart,
    "analytic_check_passed" -> passed,
    "comparison" -> "Stored returned numbers padded to 80 digits only for residual computation; no added accuracy is claimed.",
    "normalization" -> "d^D l/(i Pi^(D/2)) per loop; no exp(EulerGamma eps) factor",
    "scope" -> "Original automatic workflow smoke test against an analytic formula; not a matched Rust timing."
  |>, "RawJSON"];
  assert[passed, "Analytic agreement is below 20 digits."];
  assert[Kernels[] === {}, "Unexpected parallel kernels remain."];
  Print["ORACLE PASS: ", oracleCase, " ", mode, "; max scaled error = ",
    ToString[Max[errors], InputForm]];
  Exit[0],
  fail["Original AMFlow aborted."]
];
fail["Oracle reached the end without an explicit successful exit."];
