(* Mathematical data extraction only. No upstream solver implementation is loaded. *)
SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> 1];
args = If[Length[$ScriptCommandLine] > 0, Rest[$ScriptCommandLine],
  scriptPosition = FirstPosition[$CommandLine, "-script", Missing["none"], {1}];
  If[MissingQ[scriptPosition], {}, Take[Drop[$CommandLine, First[scriptPosition]+1], UpTo[2]]]];
If[Length[args] != 2, Print["Expected canonical_basis_NP_system.m and output JSON"]; Exit[1]];
{input, outputPath} = args;
If[FileHash[input,"SHA256","HexString"] != "cf673b527f1fcd25fc2df9f6c04998d5027b0909cd9d64b32e758d86eabbf39c", Print["Wrong pinned primary input"]; Exit[1]];
v = Get[input];
If[!ListQ[v] || Length[v] != 64, Print["unexpected source dimensions"]; Exit[1]];
norm = (-sh)^(2 eps);
plain = v /. norm -> 1;
ints = DeleteDuplicates[Cases[plain, (PL | PLx12 | PLx123 | NP)[_List, _], Infinity]];
If[!And @@ (Last[#] === 4 - 2 eps & /@ ints), Print["unexpected integral dimension"]; Exit[2]];
roots = DeleteDuplicates[Cases[plain, Power[a_, q_Rational] /; Denominator[q] == 2 :> a, Infinity]];
rootSymbols = Table[Symbol["root" <> ToString[i - 1]], {i, Length[roots]}];
rootRules = Power[a_, q_Rational] /; Denominator[q] == 2 :> rootSymbols[[First@FirstPosition[roots,a,Missing["none"],{1},Heads->False]]]^(2q);
weights = Table[Together[Coefficient[Expand[plain[[i]]], ints[[j]]] /. rootRules],
  {i, Length[plain]}, {j, Length[ints]}];
If[!And @@ Table[Together[(plain[[i]] /. rootRules) - weights[[i]].ints] === 0,
    {i, Length[plain]}], Print["linear reconstruction failed"]; Exit[3]];
rootReconstruction = And @@ Table[
  Together[((weights[[i]] /. Thread[rootSymbols -> (Sqrt /@ roots)]).ints) - plain[[i]]] === 0,
  {i, Length[plain]}];
If[!rootReconstruction, Print["root reconstruction failed"]; Exit[4]];
Print["exact original-root roundtrip=", rootReconstruction];
str[x_] := StringReplace[ToString[x, InputForm], {"\\[Rho]" -> "rho", "ρ" -> "rho", " " -> ""}];
data = <|"schema" -> "rustflow-published-gg-hg-canonical-v1",
  "paper" -> "2007.09813v2", "dimension" -> Length[v],
  "normalization" -> "(-mh2)^(2*eps)/Gamma(1+eps)^2",
  "variables" -> <|"sh" -> "mh2", "y" -> "-t/mh2", "z" -> "-u/mh2", "rho" -> "-mv2/mh2"|>,
  "roots" -> MapThread[<|"name" -> ToString[#1, InputForm], "radicand" -> str[#2]|> &, {rootSymbols, roots}],
  "integrals" -> (<|"head" -> ToString[Head[#]], "powers" -> First[#], "dimension" -> 4|> & /@ ints),
  "rows" -> Table[Select[Table[<|"integral" -> j - 1, "coefficient" -> str[weights[[i,j]]]|>,
      {j, Length[ints]}], #["coefficient"] != "0" &], {i, Length[v]}],
  "extraction" -> <|"exact_linear_reconstruction" -> True, "exact_root_reconstruction" -> rootReconstruction, "source_sha256" -> FileHash[input, "SHA256", "HexString"]|>|>;
data["provenance"] = <|"url" -> "https://arxiv.org/src/2007.09813v2/anc/differential_equations/canonical_basis_NP_system.m", "sha256" -> FileHash[input,"SHA256","HexString"], "interpretation" -> "Published 64-component basis only; no equivalence with plugin 48/61 basis established. PLx12 and PLx123 crossing definitions are not inferred from their labels."|>;
Export[outputPath, data, "RawJSON"];
Print["rows=", Length[v], " integrals=", Length[ints], " roots=", Length[roots], " entries=", Count[weights, x_ /; x =!= 0, {2}]];
Exit[0];
