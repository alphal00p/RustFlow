(* Standalone oracle for the G[1,0,1;4] example in DiffExp's
   MultiplePolylogarithms.nb. This independently written driver loads the
   upstream package unchanged. It does not supply a boundary to the Rust test.

   DIFFEXP_ROOT=/path/to/diffexp DIFFEXP_OUTPUT=/new/output/directory \
   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
   timeout --signal=TERM --kill-after=10s 120s \
   bern-wolfram -noinit -noprompt -script scripts/diffexp_mpl_oracle.wl

   DIFFEXP_OUTPUT must not exist. The source hash pins the tested package;
   the checkout is never modified. Timings cover TransportTo only and are
   single-run diagnostics, not a performance comparison. *)

$HistoryLength = 0;
driverFile = ExpandFileName[$InputFileName];
fail[msg_] := (Print["ORACLE FAILURE: ", msg]; Exit[1]);
requiredEnvironment[name_] := Module[{value = Environment[name]},
  If[!StringQ[value] || StringLength[value] == 0,
    fail["Set " <> name <> " before running this script"]];
  ExpandFileName[value]
];
sha256[path_] := IntegerString[FileHash[path, "SHA256"], 16, 64];
upstreamRoot = requiredEnvironment["DIFFEXP_ROOT"];
outputRoot = requiredEnvironment["DIFFEXP_OUTPUT"];
package = FileNameJoin[{upstreamRoot, "DiffExp.m"}];
expectedHash = "67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7";
If[!FileExistsQ[package], fail["DIFFEXP_ROOT does not contain DiffExp.m"]];
If[sha256[package] =!= expectedHash, fail["Pinned DiffExp.m hash mismatch"]];
If[FileExistsQ[outputRoot] || DirectoryQ[outputRoot],
  fail["DIFFEXP_OUTPUT must be a new directory"]];
If[CreateDirectory[FileNameJoin[{outputRoot, "matrices"}],
    CreateIntermediateDirectories -> True] === $Failed,
  fail["Could not create output directory"]];
Get[package];

(* Basis: {G[1,0,1;t], G[0,1;t], G[1;t], 1}. All inputs stay exact. *)
matrix = {{0, 1/(t-1), 0, 0}, {0, 0, 1/t, 0},
          {0, 0, 0, 1/(t-1)}, {0, 0, 0, 0}};
matrixFile = FileNameJoin[{outputRoot, "matrices", "dt_0.m"}];
Put[matrix, matrixFile];
records = {};
results = {};
scalarRecord[value_] := <|
  "real_wolfram_input" -> ToString[Re[value], InputForm],
  "imaginary_wolfram_input" -> ToString[Im[value], InputForm],
  "precision_wolfram_input" -> ToString[Precision[value], InputForm]|>;

CheckAbort[
  Do[
    DiffExp`LoadConfiguration[{
      DiffExp`EpsilonOrder -> 0, WorkingPrecision -> 250,
      DiffExp`ChopPrecision -> 225, DiffExp`ExpansionOrder -> order,
      DiffExp`UseMobius -> True, DiffExp`UsePade -> False,
      DiffExp`DivisionOrder -> 3, DiffExp`Verbosity -> 1,
      DiffExp`MatrixDirectory -> FileNameJoin[{outputRoot, "matrices"}],
      DiffExp`DeltaPrescriptions -> {t-1-I*\[Delta], t+I*\[Delta]},
      "Parallel" -> False}];
    boundary = DiffExp`PrepareBoundaryConditions[
      {0+O[x]^(1/2), 0+O[x]^(1/2), 0+O[x]^(1/2), 1+O[x]^(1/2)},
      {t -> 4*x}];
    measured = AbsoluteTiming[DiffExp`TransportTo[boundary, {t -> 4}]];
    result = measured[[2]];
    If[!ListQ[result] || Length[result] < 3 ||
       Dimensions[result[[2]]] =!= {4, 1} ||
       !AllTrue[Flatten[result[[2]]], NumberQ],
      fail["Malformed returned vector"]];
    analyticErrors = Abs[result[[2, 2;;4, 1]] -
      N[{-PolyLog[2, 4], Log[3]+I*Pi, 1}, 250]];
    If[!TrueQ[Max[analyticErrors] < 10^-20],
      fail[{"Analytic subchain disagreement", analyticErrors}]];
    rawFile = FileNameJoin[{outputRoot, "result-order" <> ToString[order] <> ".wl"}];
    Put[result, rawFile];
    AppendTo[results, result];
    AppendTo[records, <|
      "expansion_order" -> order,
      "transport_wall_seconds" -> measured[[1]],
      "values" -> (scalarRecord /@ result[[2, All, 1]]),
      "reported_errors_wolfram_input" -> (ToString[#, InputForm]& /@ Flatten[result[[3]]]),
      "analytic_subchain_errors_wolfram_input" -> (ToString[#, InputForm]& /@ analyticErrors),
      "raw_result_file" -> FileNameTake[rawFile],
      "raw_result_sha256" -> sha256[rawFile]|>],
    {order, {50, 75}}];

  refinementErrors = Abs[results[[1, 2, All, 1]] - results[[2, 2, All, 1]]];
  notebookReference =
    -6.7782180257804207212554826775005988168291802221955692129682`36.14650810661106 +
    I*0.9250147943833369547396749852220309435917997631163983727603`35.281541251279386;
  referenceError = Abs[results[[2, 2, 1, 1]] - notebookReference];
  If[!TrueQ[Max[refinementErrors] < 10^-20 && referenceError < 10^-30],
    fail[{"Refinement/reference disagreement", refinementErrors, referenceError}]];
  metadata = <|
    "status" -> "passed",
    "upstream_commit" -> "784c8229bf92369a03f011a48e161522c8c54bbd",
    "source_sha256" -> sha256[package], "driver_sha256" -> sha256[driverFile],
    "matrix_sha256" -> sha256[matrixFile], "wolfram_version" -> $Version,
    "case" -> "G[1,0,1;4]", "basis" -> {"G[1,0,1;t]", "G[0,1;t]", "G[1;t]", "1"},
    "exact_boundary" -> {"t -> 0+", {"0", "0", "0", "1"}},
    "exact_endpoint" -> "t -> 4", "initial_parameterization" -> "t -> 4*x",
    "prescriptions" -> {"t-1-I*delta", "t+I*delta"},
    "settings" -> <|"working_precision" -> 250, "chop_precision" -> 225,
      "epsilon_order" -> 0, "expansion_orders" -> {50, 75},
      "division_order" -> 3, "mobius" -> True, "pade" -> False,
      "parallel" -> False, "verbosity" -> 1|>,
    "notebook_deviations" -> {"Verbosity is 1 to preserve diagnostics (notebook 0).",
      "Duplicate identical t-1 delta prescription is represented once."},
    "runs" -> records,
    "absolute_order_refinement_errors_wolfram_input" -> (ToString[#, InputForm]& /@ refinementErrors),
    "recorded_notebook_reference_wolfram_input" -> ToString[notebookReference, InputForm],
    "recorded_reference_error_wolfram_input" -> ToString[referenceError, InputForm],
    "asserted_absolute_refinement_digits" -> 20,
    "asserted_absolute_reference_digits" -> 30,
    "precision_note" -> "Working precision is not an achieved-accuracy claim. Notebook reference precision is preserved.",
    "timing_scope" -> "One TransportTo timing per order; excludes launch, configuration and boundary preparation. No repeated benchmark or Rust parity claim."|>;
  If[Export[FileNameJoin[{outputRoot, "result.json"}], metadata, "RawJSON"] === $Failed,
    fail["Could not export metadata"]];
  Print["ORACLE PASS"]; Exit[0],
  fail["DiffExp aborted"]
];
